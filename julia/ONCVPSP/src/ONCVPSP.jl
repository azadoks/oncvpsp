"""
ONCVPSP's runs and stages in Julia, over liboncvpsp (include/oncvpsp.h).

The input is a dictionary in the layout of the TOML input, the TOML text,
or the path of a TOML file; keyword arguments are tables merged over it,
key by key, and nothing values are left out (TOML has no null). Each
function returns a `Handle`, whose fields are arrays in the Fortran shapes
(Julia's order is Fortran's); a failure throws an `OncvpspError`, or with
`check=false` is left in the handle. `text=true` keeps the text output
(`h.text`), `false` discards it, a path writes it there. The library is `ENV["ONCVPSP_LIBRARY"]`, or liboncvpsp on the
system's search path.

    using ONCVPSP
    run = ONCVPSP.run("tests/data/toml/14_Si.toml"; irel=1)
    size(run["pot.vkb"])   # (mmax, mxprj, lmax+1, nkap)
    ONCVPSP.write(run, "Si.upf")
"""
module ONCVPSP

using Libdl
using TOML

export Handle, OncvpspError, all_electron, build_info, channel, confine, coupled, has_node

const API_VERSION = 4
const MXLIN = 5
const MXRANK = 5
const MSG_PSPOT_NODE = 122
const FAILURES = ("none", "input", "convergence", "not_found", "construction",
                  "numerical", "io", "internal")
const STAGES = ("none", "input", "ae", "channel", "vkb", "valence", "model_core",
                "screening", "diag", "ghosts", "unscreen", "configs", "plot",
                "phsft", "tests", "output", "confine")

"oncv_failure_t"
struct Failure
    kind::Cint
    msg::Cint
    ierr::Cint
    stage::Cint
    ll::Cint
    ikap::Cint
    iproj::Cint
end

"oncv_dims_t"
struct Dims
    version::Cint
    irel::Cint
    nkap::Cint
    iexc::Cint
    mmax::Cint
    lmax::Cint
    lloc::Cint
    nproj::NTuple{MXLIN + 1,Cint}
    nrl::Cint
end

const libhandle = Ref{Ptr{Cvoid}}(C_NULL)
const symbols = Dict{Symbol,Ptr{Cvoid}}()
const fields = Dict{String,Char}()

"the address of the library's function name"
sym(name::Symbol) = get!(() -> Libdl.dlsym(libhandle[], name), symbols, name)

function __init__()
    path = get(ENV, "ONCVPSP_LIBRARY", "")
    if isempty(path)
        path = Libdl.find_library(["liboncvpsp", "oncvpsp"])
    end
    isempty(path) && error("liboncvpsp not found; set ONCVPSP_LIBRARY to its path")
    # global, as MKL's dispatch libraries need the symbols of the MKL it links
    libhandle[] = Libdl.dlopen(path, Libdl.RTLD_GLOBAL | Libdl.RTLD_NOW)
    empty!(symbols)
    version = ccall(sym(:oncv_api_version), Cint, ())
    version == API_VERSION ||
        error("liboncvpsp $path has interface version $version, this package $API_VERSION")
    buf = Vector{UInt8}(undef, 64)
    empty!(fields)
    while true
        nn = ccall(sym(:oncv_field_name), Cint, (Cint, Ptr{UInt8}, Csize_t),
                   length(fields) + 1, buf, length(buf))
        nn == 0 && break
        kind, name = split(String(buf[1:nn]))
        fields[name] = kind[1]
    end
end

"""
the library's build: a Dict of "commit", "source" and "compiled", its git
commit, the state of the source against it and the build time
"""
function build_info()
    nout = Ref{Csize_t}(0)
    ccall(sym(:oncv_build_info), Cint, (Ptr{UInt8}, Csize_t, Ref{Csize_t}), C_NULL, 0, nout)
    buf = Vector{UInt8}(undef, nout[])
    ccall(sym(:oncv_build_info), Cint, (Ptr{UInt8}, Csize_t, Ref{Csize_t}),
          buf, length(buf), nout)
    lines = filter(!isempty, split(String(buf), '\n'))
    Dict(String(kk) => String(vv) for (kk, vv) in (split(ll, ' '; limit=2) for ll in lines))
end

"""
A failed run or stage: ifail and its parts (oncv_failure_t), the names of
its kind and stage, and the text output (the message is there) if kept
"""
struct OncvpspError <: Exception
    ifail::Int
    failure::Failure
    text::Union{String,Nothing}
end

kind_name(err::OncvpspError) = get(FAILURES, err.failure.kind + 1, "?")
stage_name(err::OncvpspError) = get(STAGES, err.failure.stage + 1, "?")

function Base.showerror(io::IO, err::OncvpspError)
    ff = err.failure
    print(io, "$(kind_name(err)) failure (ifail $(err.ifail), msg $(ff.msg)) at stage ",
          stage_name(err), ff.ll >= 0 ? ", l=$(ff.ll)" : "")
    errors = filter(ll -> occursin("ERROR", ll), split(something(err.text, ""), '\n'))
    isempty(errors) || print(io, ": ", strip(errors[end]))
end

# the fields of each kind of handle, by prefix ("" for the run's); the
# library writes an error for each field asked of a handle that holds none
const STAGE_PREFIXES = ("ae.", "channel.", "confine.")
const HOLDS = Dict(:run => ("", "ae."), :ae => ("ae.",), :channel => ("channel.",),
                   :confine => ("confine.",), :none => ())

function prefix(name)
    ii = findfirst(pp -> startswith(name, pp), STAGE_PREFIXES)
    ii === nothing ? "" : STAGE_PREFIXES[ii]
end

"""
One calculation of liboncvpsp, of a kind (:run, :ae, :channel, :confine):
`h[name]` is a field, `keys(h)` those it holds. The handle is freed by a
finalizer, or by `close(h)`
"""
mutable struct Handle
    ptr::Ptr{Cvoid}
    kind::Symbol
    holds::Tuple
    keep_text::Bool
    ifail::Int
    function Handle(kind::Symbol, text=true)
        ptr = ccall(sym(:oncv_new), Ptr{Cvoid}, ())
        ptr == C_NULL && throw(OutOfMemoryError())
        hh = new(ptr, kind, HOLDS[kind], text === true, 0)
        finalizer(close, hh)
        if text === true
            ccall(sym(:oncv_set_output_buffer), Cint, (Ptr{Cvoid}, Cint), ptr, 1)
        elseif text !== false
            path = String(text)
            ccall(sym(:oncv_set_output), Cint, (Ptr{Cvoid}, Ptr{UInt8}, Csize_t),
                  ptr, path, sizeof(path))
        end
        hh
    end
end

"frees the handle now"
function Base.close(hh::Handle)
    if hh.ptr != C_NULL
        ccall(sym(:oncv_free), Cvoid, (Ptr{Cvoid},), hh.ptr)
        hh.ptr = C_NULL
    end
    nothing
end

function ptr(hh::Handle)
    hh.ptr == C_NULL && throw(ArgumentError("the handle is closed"))
    hh.ptr
end

function done!(hh::Handle, ifail, check)
    hh.ifail = ifail
    ifail != 0 && (hh.holds = ())
    ifail != 0 && check && throw(error_of(hh))
    hh
end

"the OncvpspError of the handle's last run or stage"
error_of(hh::Handle) = OncvpspError(hh.ifail, failure(hh), text(hh))

"the failure record (oncv_failure_t) of the last run or stage"
function failure(hh::Handle)
    ff = Ref{Failure}()
    ccall(sym(:oncv_failure), Cint, (Ptr{Cvoid}, Ref{Failure}), ptr(hh), ff)
    ff[]
end

"the dimensions (oncv_dims_t), checked of this interface's version"
function dims(hh::Handle)
    dd = Ref{Dims}()
    ccall(sym(:oncv_dims), Cint, (Ptr{Cvoid}, Ref{Dims}), ptr(hh), dd)
    dd[].version == API_VERSION || error("oncv_dims has version $(dd[].version)")
    dd[]
end

"the text output of the last run or stage, nothing if not kept"
function text(hh::Handle)
    hh.keep_text || return nothing
    nout = Ref{Csize_t}(0)
    ccall(sym(:oncv_output_text), Cint, (Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Ref{Csize_t}),
          ptr(hh), C_NULL, 0, nout)
    buf = Vector{UInt8}(undef, nout[])
    ccall(sym(:oncv_output_text), Cint, (Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Ref{Csize_t}),
          ptr(hh), buf, length(buf), nout)
    String(buf)
end

function Base.getproperty(hh::Handle, name::Symbol)
    name === :text && return text(hh)
    name === :dims && return dims(hh)
    name === :failure && return failure(hh)
    getfield(hh, name)
end

"the shape of a field (Fortran order), nothing if the handle holds none"
function shape(hh::Handle, name::AbstractString)
    (haskey(fields, name) && prefix(name) in hh.holds) || return nothing
    rank, shp = Ref{Cint}(0), zeros(Cint, MXRANK)
    ifail = ccall(sym(:oncv_shape), Cint,
                  (Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Ref{Cint}, Ptr{Cint}),
                  ptr(hh), name, sizeof(name), rank, shp)
    ifail == 0 ? Tuple(Int.(shp[1:rank[]])) : nothing
end

function Base.getindex(hh::Handle, name::AbstractString)
    shp = shape(hh, name)
    shp === nothing && throw(KeyError(name))
    if fields[name] == 'r'
        out = Array{Float64}(undef, shp...)
        ifail = ccall(sym(:oncv_get_double), Cint,
                      (Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Ptr{Float64}, Csize_t),
                      ptr(hh), name, sizeof(name), out, length(out))
    else
        out = Array{Cint}(undef, shp...)
        ifail = ccall(sym(:oncv_get_int), Cint,
                      (Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Ptr{Cint}, Csize_t),
                      ptr(hh), name, sizeof(name), out, length(out))
    end
    ifail == 0 || error("oncv_get of $name: ifail $ifail")
    shp == () ? out[] : out
end

Base.haskey(hh::Handle, name::AbstractString) = shape(hh, name) !== nothing
Base.keys(hh::Handle) = sort!([nn for nn in keys(fields) if haskey(hh, nn)])

"""
writes the pseudopotential file of a run: psp8, upf or psml, by default
from path's extension
"""
function write(hh::Handle, path::AbstractString; format=lstrip(splitext(path)[2], '.'))
    fmt = String(format)
    ifail = ccall(sym(:oncv_write_files), Cint,
                  (Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t),
                  ptr(hh), fmt, sizeof(fmt), path, sizeof(path))
    ifail == 0 || (hh.ifail = ifail; throw(error_of(hh)))
    nothing
end

Base.show(io::IO, hh::Handle) =
    print(io, "ONCVPSP.Handle(:$(hh.kind), ",
          hh.ptr == C_NULL ? "closed" : "ifail=$(hh.ifail)", ")")

# the input

"""
the TOML text of the input `inp` (a dictionary in the TOML layout, the
TOML text, or the path of a TOML file: a string of one line ending in
.toml) with the tables merged over it, key by key
"""
function input_text(inp, tables)
    if inp isa AbstractString && !occursin('\n', inp) && endswith(inp, ".toml")
        inp = read(inp, String)
    end
    if inp isa AbstractString
        isempty(tables) && return String(inp)
        inp = TOML.parse(inp)
    end
    data = Dict{String,Any}(String(kk) => (vv isa AbstractDict ? Dict{String,Any}(vv) : vv)
                            for (kk, vv) in something(inp, Dict()))
    for (kk, vv) in tables
        key = String(kk)
        if vv isa AbstractDict && get(data, key, nothing) isa AbstractDict
            merge!(data[key], Dict{String,Any}(String(k) => v for (k, v) in vv))
        else
            data[key] = vv
        end
    end
    sprint(TOML.print, dropnothing(data))
end

"the input with its nothing values left out (TOML has no null)"
dropnothing(dd::AbstractDict) =
    Dict{String,Any}(String(kk) => dropnothing(vv) for (kk, vv) in dd if vv !== nothing)
dropnothing(vv::AbstractVector) = [dropnothing(xx) for xx in vv]
dropnothing(vv) = vv

irel_of(irel::Integer) = irel in 0:2 ? Cint(irel) : throw(ArgumentError("irel is 0, 1 or 2"))
irel_of(irel::Union{Symbol,AbstractString}) =
    Cint(Dict("nr" => 0, "sr" => 1, "fr" => 2, "r" => 2)[String(irel)])

"""
the whole run of the input (oncv_run_toml), or from the all-electron atom
of `ae` (all_electron, or a run), bit for bit the same (oncv_run_ae_toml);
`skip_ld` leaves out the log derivatives
"""
function run(inp=nothing; irel=1, ae=nothing, skip_ld=false, check=true, text=true,
             tables...)
    data = input_text(inp, tables)
    hh = Handle(:run, text)
    ifail = if ae === nothing
        ccall(sym(:oncv_run_toml), Cint, (Ptr{Cvoid}, Cint, Ptr{UInt8}, Csize_t, Cint),
              ptr(hh), irel_of(irel), data, sizeof(data), skip_ld)
    else
        ccall(sym(:oncv_run_ae_toml), Cint,
              (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Cint),
              ptr(hh), ptr(ae), data, sizeof(data), skip_ld)
    end
    done!(hh, ifail, check)
end

"the all-electron atom of the input (oncv_ae_toml): fields ae.*"
function all_electron(inp=nothing; irel=1, check=true, text=true, tables...)
    data = input_text(inp, tables)
    hh = Handle(:ae, text)
    ifail = ccall(sym(:oncv_ae_toml), Cint, (Ptr{Cvoid}, Cint, Ptr{UInt8}, Csize_t),
                  ptr(hh), irel_of(irel), data, sizeof(data))
    done!(hh, ifail, check)
end

const MODES = Dict(:channel => 0, :projectors => 1, :node => 2)

"""
the channel of angular momentum `l` of the input, from the all-electron
atom of `ae` (oncv_channel_toml): fields channel.*. `mode` :projectors
also makes the local potential and the channel's projectors; :node stops
at the node test of the first projector (has_node)
"""
function channel(inp, ae::Handle, l::Integer; mode=:channel, check=true, text=true,
                 tables...)
    data = input_text(inp, tables)
    hh = Handle(mode === :node ? :none : :channel, text)
    ifail = ccall(sym(:oncv_channel_toml), Cint,
                  (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Cint, Cint),
                  ptr(hh), ptr(ae), data, sizeof(data), l, MODES[mode])
    done!(hh, ifail, check)
end

"""
whether the first pseudo-wave of channel `l` has a node (pspot's test),
from the all-electron atom of `ae`; a failure before the test throws
"""
function has_node(inp, ae::Handle, l::Integer; text=false, tables...)
    hh = channel(inp, ae, l; mode=:node, check=false, text=text, tables...)
    hh.ifail == 0 && return false
    failure(hh).msg == MSG_PSPOT_NODE && return true
    throw(error_of(hh))
end

"""
the model core, screening and tests again on the stages of `base`'s run
with the model core icmod, fcfact, rcfact (oncv_coupled): a run
"""
function coupled(base::Handle, icmod::Integer, fcfact=0.0, rcfact=0.0; skip_ld=false,
                 check=true, text=true)
    hh = Handle(:run, text)
    ifail = ccall(sym(:oncv_coupled), Cint,
                  (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Cdouble, Cdouble, Cint),
                  ptr(hh), ptr(base), icmod, fcfact, rcfact, skip_ld)
    done!(hh, ifail, check)
end

"""
the confined atom of `base`'s run, as oncvpsp_confine.x (oncv_confine):
fields confine.*. `radii` of the sphere (bohr), nothing for the default
(1.25, 1.5, 2 times the largest core radius), () for no sphere;
`channels`' l, nothing for all; `energies` (emin, emax, de) in Ha, nothing
for the input's window and 0.1; `walls` (Ha/bohr**2), nothing for the
default, () for none; `r0` the walls' onset, `ecutmax` the basis's largest
cutoff (Ha), 0 for their defaults
"""
function confine(base::Handle; radii=nothing, channels=nothing, energies=nothing,
                 walls=nothing, r0=0.0, ecutmax=0.0, check=true, text=true)
    rr = radii === nothing ? Float64[] : Float64[radii...]
    nradii = radii === nothing ? 0 : (isempty(rr) ? -1 : length(rr))
    cc = channels === nothing ? Cint[] : Cint[channels...]
    energies === nothing || length(energies) == 3 ||
        throw(ArgumentError("energies is (emin, emax, de)"))
    ee = energies === nothing ? C_NULL : Float64[energies...]
    ww = walls === nothing ? Float64[] : Float64[walls...]
    nwalls = walls === nothing ? -1 : length(ww)
    hh = Handle(:confine, text)
    ifail = GC.@preserve rr cc ee ww ccall(sym(:oncv_confine), Cint,
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Float64}, Cint, Ptr{Cint}, Cint, Ptr{Float64},
         Ptr{Float64}, Cint, Cdouble, Cdouble),
        ptr(hh), ptr(base), rr, nradii, cc, length(cc), ee, ww, nwalls, r0, ecutmax)
    done!(hh, ifail, check)
end

end # module
