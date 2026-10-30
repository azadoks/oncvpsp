"""
ONCVPSP's runs and stages in Julia, over liboncvpsp (include/oncvpsp.h).

The input is a dictionary in the layout of the TOML input, the TOML text,
or the path of a TOML file; keyword arguments are tables merged over it,
key by key, and nothing values are left out (TOML has no null). Each
function returns a handle (`Run`, `AllElectron`, `ONCVPSP.Channel`,
`Confinement`), an `AbstractDict` of its fields to arrays in the Fortran
shapes (Julia's order is Fortran's), which `fieldinfo` describes; a
failure throws an `OncvpspError`, or with `check=false` is left in the
handle. `text=true` keeps the text output (`h.text`), `false` discards
it, a path writes it there. The library is `ENV["ONCVPSP_LIBRARY"]`, or
liboncvpsp on the system's search path.

    using ONCVPSP
    oncvpsp("tests/data/toml/14_Si.toml"; irel=:sr) do run
        size(run["pot.vkb"])   # (mmax, mxprj, mxl+1, nkap)
        kernel(run, 0)         # the nonlocal operator of l=0, (mmax, mmax)
        write("Si.upf", run)
    end
"""
module ONCVPSP

using Libdl
using TOML

export AllElectron, Confinement, Handle, Local, OncvpspError, Run, all_electron, assemble,
       build_info, channel, confine, convergence, core_density, coupled, fieldinfo, has_node,
       kernel, local_potential, mesh, model_core, oncvpsp, projectors, semilocal,
       valence_density

const API_VERSION = 9
const MXLIN = 5
const MXRANK = 5
const MSG_PSPOT_NODE = 122
const FAILURES = (:none, :input, :convergence, :not_found, :construction, :numerical,
                  :io, :internal)
const STAGES = (:none, :input, :ae, :channel, :vkb, :valence, :model_core, :screening,
                :diag, :ghosts, :unscreen, :configs, :plot, :phsft, :tests, :output,
                :confine)

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

"""
A field of the results (oncv_field_info): its name, kind ('r' real, 'i'
integer), axes, the first index first (an axis of different quantities as
their names in brackets, "[ecut,err]"), units ("" for none) and what it is
"""
struct FieldInfo
    name::String
    kind::Char
    axes::Vector{String}
    units::String
    doc::String
end

function Base.show(io::IO, ::MIME"text/plain", ff::FieldInfo)
    axes = isempty(ff.axes) ? "scalar" : join(ff.axes, ", ")
    print(io, ff.name, " (", axes, ")", isempty(ff.units) ? "" : " [$(ff.units)]", ": ",
          ff.doc)
end

const fields = Dict{String,FieldInfo}()
const fieldorder = String[]

"what the field `name` is: its axes, units and a sentence (a `FieldInfo`)"
fieldinfo(name::AbstractString) =
    get(() -> throw(KeyError(name)), fields, name)

function info_text(ifield, what)
    nout = Ref{Csize_t}(0)
    ccall(sym(:oncv_field_info), Cint, (Cint, Cint, Ptr{UInt8}, Csize_t, Ref{Csize_t}),
          ifield, what, C_NULL, 0, nout)
    buf = Vector{UInt8}(undef, nout[])
    ccall(sym(:oncv_field_info), Cint, (Cint, Cint, Ptr{UInt8}, Csize_t, Ref{Csize_t}),
          ifield, what, buf, length(buf), nout)
    String(buf)
end

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
    empty!(fieldorder)
    while true
        ifield = length(fieldorder) + 1
        nn = ccall(sym(:oncv_field_name), Cint, (Cint, Ptr{UInt8}, Csize_t),
                   ifield, buf, length(buf))
        nn == 0 && break
        kind, name = String.(split(String(buf[1:nn])))
        fields[name] = FieldInfo(name, kind[1], String.(split(info_text(ifield, 1))),
                                 info_text(ifield, 2), info_text(ifield, 3))
        push!(fieldorder, name)
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
A failed run or stage: ifail, its record (oncv_failure_t) and the text
output (the message is there) if kept. `err.kind` and `err.stage` are
Symbols (`:input`, `:convergence`, ...; `:ae`, `:channel`, ...), and
`err.msg`, `err.ll`, `err.ikap`, `err.iproj` the record's; `isnode(err)`
whether it is a node in a channel's first pseudo-wave (has_node)
"""
struct OncvpspError <: Exception
    ifail::Int
    failure::Failure
    text::Union{String,Nothing}
end

function Base.getproperty(err::OncvpspError, name::Symbol)
    ff = getfield(err, :failure)
    name === :kind && return get(FAILURES, ff.kind + 1, :unknown)
    name === :stage && return get(STAGES, ff.stage + 1, :unknown)
    name in (:msg, :ierr, :ll, :ikap, :iproj) && return getfield(ff, name)
    getfield(err, name)
end

Base.propertynames(::OncvpspError) =
    (:ifail, :failure, :text, :kind, :stage, :msg, :ierr, :ll, :ikap, :iproj)

isnode(err::OncvpspError) = err.msg == MSG_PSPOT_NODE
kind_name(err::OncvpspError) = String(err.kind)
stage_name(err::OncvpspError) = String(err.stage)

function Base.showerror(io::IO, err::OncvpspError)
    print(io, "$(err.kind) failure (ifail $(err.ifail), msg $(err.msg)) at stage ",
          err.stage, err.ll >= 0 ? ", l=$(err.ll)" : "")
    errors = filter(ll -> occursin("ERROR", ll), split(something(err.text, ""), '\n'))
    isempty(errors) || print(io, ": ", strip(errors[end]))
end

# the fields of each kind of handle, by prefix ("" for the run's); the
# library writes an error for each field asked of a handle that holds none
const STAGE_PREFIXES = ("ae.", "channel.", "confine.", "local.")
holds_of(K) = Dict(:run => ("", "ae."), :ae => ("ae.",), :channel => ("channel.",),
                   :confine => ("confine.",), :local => ("local.",))[K]

function prefix(name)
    ii = findfirst(pp -> startswith(name, pp), STAGE_PREFIXES)
    ii === nothing ? "" : STAGE_PREFIXES[ii]
end

"""
One calculation of liboncvpsp, of a kind K: `Run` (`Handle{:run}`),
`AllElectron`, `ONCVPSP.Channel` or `Confinement`. An `AbstractDict` of
the fields it holds: `h[name]` is a field, `keys(h)` those it holds;
`h.text`, `h.dims`, `h.failure` and `h.ifail` its run's. The handle is
freed by a finalizer, or by `close(h)`
"""
mutable struct Handle{K} <: AbstractDict{String,Any}
    ptr::Ptr{Cvoid}
    holds::Tuple
    keep_text::Bool
    ifail::Int
    function Handle{K}(text=true) where {K}
        ptr = ccall(sym(:oncv_new), Ptr{Cvoid}, ())
        ptr == C_NULL && throw(OutOfMemoryError())
        hh = new{K}(ptr, holds_of(K), text === true, 0)
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

"a whole run: pot.*, tests.*, core.*, warnings.* and its all-electron atom ae.*"
const Run = Handle{:run}
"the all-electron atom alone: ae.*"
const AllElectron = Handle{:ae}
"one channel: channel.* (not exported: Base has a Channel)"
const Channel = Handle{:channel}
"the confined atom of a run: confine.*"
const Confinement = Handle{:confine}
"the local potential alone: local.*, from which channels and runs are assembled"
const Local = Handle{:local}

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

Base.propertynames(::Handle) = (:text, :dims, :failure, :ifail)

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
    if fields[name].kind == 'r'
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
Base.keys(hh::Handle) = [nn for nn in fieldorder if haskey(hh, nn)]
Base.length(hh::Handle) = count(nn -> haskey(hh, nn), fieldorder)
Base.get(hh::Handle, name::AbstractString, default) = haskey(hh, name) ? hh[name] : default
function Base.iterate(hh::Handle, state=(keys(hh), 1))
    names, ii = state
    ii > length(names) ? nothing : (names[ii] => hh[names[ii]], (names, ii + 1))
end
# one calculation: equal only to itself (not AbstractDict's equal fields)
Base.:(==)(aa::Handle, bb::Handle) = aa === bb
Base.hash(hh::Handle, hh0::UInt) = hash(objectid(hh), hh0)

"""
writes the pseudopotential file of a run: psp8, upf or psml, by default
from path's extension
"""
function Base.write(path::AbstractString, hh::Run; format=lstrip(splitext(path)[2], '.'))
    fmt = String(format)
    ifail = ccall(sym(:oncv_write_files), Cint,
                  (Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t),
                  ptr(hh), fmt, sizeof(fmt), path, sizeof(path))
    ifail == 0 || (hh.ifail = ifail; throw(error_of(hh)))
    nothing
end

function Base.show(io::IO, hh::Handle{K}) where {K}
    name = Dict(:run => "Run", :ae => "AllElectron", :channel => "ONCVPSP.Channel",
                :confine => "Confinement", :local => "Local")[K]
    state = hh.ptr == C_NULL ? "closed" :
            hh.ifail != 0 ? "ifail=$(hh.ifail)" : "$(length(hh)) fields"
    print(io, name, "(", state, ")")
end
Base.show(io::IO, ::MIME"text/plain", hh::Handle) = show(io, hh)

# the physics of the fields

function kap(ikap, nkap)
    1 <= ikap <= nkap || throw(ArgumentError("ikap is 1 to $nkap, not $ikap"))
    ikap
end

function lcheck(l, lmax)
    0 <= l <= lmax || throw(ArgumentError("l is 0 to $lmax, not $l"))
    l
end

"the log radial mesh (bohr)"
mesh(run::Run) = run["pot.rr"]
mesh(ae::AllElectron) = ae["ae.rr"]

"the local potential (Ha) on the mesh; of a `Local`, before a run's ionic tail"
local_potential(run::Run) = run["pot.vpuns"][:, run.dims.lloc+1]
local_potential(loc::Local) = loc["local.vloc"]

"the unscreened semi-local potential of `l` (Ha) on the mesh"
semilocal(run::Run, l::Integer) = run["pot.vpuns"][:, lcheck(l, run.dims.lmax)+1]

"""
the projectors of `l` (`ikap` 2 for j=l-1/2 when fully relativistic),
nproj of them: `(energies, vkb)`, their coefficients (Ha) and the
projectors r*beta on the mesh (mmax, nproj), normalized. Only their
kernel is independent of their basis
"""
function projectors(run::Run, l::Integer; ikap::Integer=1)
    dd = run.dims
    kk, nn = kap(ikap, dd.nkap), dd.nproj[lcheck(l, dd.lmax)+1]
    (energies=run["pot.evkb"][1:nn, l+1, kk], vkb=run["pot.vkb"][:, 1:nn, l+1, kk])
end

"""
with tests, whether the channel's log derivatives are its input's whole run's
bit for bit (they read the local potential no further than rcmax, where the
run's ionic tail starts); the ghost tests always reach beyond it
"""
isexact(ch::Channel) = ch["channel.logder.exact"] != 0

"with mode :projectors, the channel's projectors, as those of a run"
function projectors(ch::Channel; ikap::Integer=1)
    kk = kap(ikap, ch["channel.mkap"])
    evkb = ch["channel.evkb"]
    isempty(evkb) &&
        throw(ArgumentError("the channel has no projectors without mode :projectors"))
    nn = ch["channel.nproj"]
    (energies=evkb[1:nn, kk], vkb=ch["channel.vkb"][:, 1:nn, kk])
end

"""
the nonlocal operator of `l` on the mesh, sum_i vkb_i(r) e_i vkb_i(r')
(Ha/bohr), (mmax, mmax): of the projectors, what does not depend on their
basis
"""
kernel(run::Run, l::Integer; ikap::Integer=1) = kernel_of(projectors(run, l; ikap=ikap))
kernel(ch::Channel; ikap::Integer=1) = kernel_of(projectors(ch; ikap=ikap))
kernel_of(pp) = pp.vkb * (pp.energies .* transpose(pp.vkb))

"""
the pseudo valence density (e/bohr^3) of the reference configuration:
∫ 4π r² ρ dr is the valence charge (the field pot.rho is 4π times it)
"""
valence_density(run::Run) = run["pot.rho"] ./ 4π

"the model core density (e/bohr^3), as `valence_density`"
model_core(run::Run) = run["pot.rhomod"][:, 1] ./ 4π

"the all-electron core density (e/bohr^3), as `valence_density`"
core_density(run::Run) = run["rhoc"] ./ 4π

"""
the convergence of the pseudo-wave `iproj` of `l` with the cutoff: rows
(Ecut, err) (Ha) of the levels it reaches, err from 1e-2 down by half
decades
"""
function convergence(run::Run, l::Integer; ikap::Integer=1, iproj::Integer=1)
    dd = run.dims
    kk, np = kap(ikap, dd.nkap), dd.nproj[lcheck(l, dd.lmax)+1]
    1 <= iproj <= np || throw(ArgumentError("iproj is 1 to $np, not $iproj"))
    table = permutedims(run["tests.convergence.cross"][:, :, iproj, l+1, kk])
    table[table[:, 1] .> 0, :]
end

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
function oncvpsp(inp=nothing; irel=1, ae=nothing, skip_ld=false, check=true, text=true,
                 tables...)
    data = input_text(inp, tables)
    hh = Run(text)
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

"`oncvpsp`, by the name of the C function and the Python package's"
const run = oncvpsp

"""
the handle of `make(args...; kwargs...)` given to `f`, closed after it
returns or throws: `oncvpsp(inp; irel=1) do run ... end`
"""
function with_handle(f, make, args...; kwargs...)
    hh = make(args...; kwargs...)
    try
        f(hh)
    finally
        close(hh)
    end
end

oncvpsp(f::Function, inp=nothing; kwargs...) = with_handle(f, oncvpsp, inp; kwargs...)
all_electron(f::Function, inp=nothing; kwargs...) =
    with_handle(f, all_electron, inp; kwargs...)

"the all-electron atom of the input (oncv_ae_toml): fields ae.*"
function all_electron(inp=nothing; irel=1, check=true, text=true, tables...)
    data = input_text(inp, tables)
    hh = AllElectron(text)
    ifail = ccall(sym(:oncv_ae_toml), Cint, (Ptr{Cvoid}, Cint, Ptr{UInt8}, Csize_t),
                  ptr(hh), irel_of(irel), data, sizeof(data))
    done!(hh, ifail, check)
end

const MODES = Dict(:channel => 0, :projectors => 1, :node => 2)

"""
the channel of angular momentum `l` of the input, from the all-electron
atom of `ae` (oncv_channel_toml): fields channel.*. `mode` :projectors
also makes the local potential (or takes `loc`'s, `local_potential`)
and the channel's projectors; :node stops at the node test of the first
projector (has_node). `tests=true`, with mode :projectors, also makes its
log derivatives (channel.logder.*, on `energies` (emin, emax, de) in Ha,
nothing for the input's window) and its ghost tests (channel.ghosts.*),
on the local potential before a whole run's ionic tail beyond rcmax:
`ONCVPSP.isexact(ch)` whether the log derivatives are the run's bit for bit
"""
function channel(inp, ae::Handle, l::Integer; mode=:channel, loc=nothing, tests=false,
                 energies=nothing, check=true, text=true, tables...)
    data = input_text(inp, tables)
    hh = Channel(text)
    mode === :node && (hh.holds = ())
    ifail = if loc === nothing && !tests
        ccall(sym(:oncv_channel_toml), Cint,
              (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Cint, Cint),
              ptr(hh), ptr(ae), data, sizeof(data), l, MODES[mode])
    else
        energies === nothing || length(energies) == 3 ||
            throw(ArgumentError("energies is (emin, emax, de)"))
        ee = energies === nothing ? C_NULL : Float64[energies...]
        GC.@preserve ee ccall(sym(:oncv_channel_local_toml), Cint,
              (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Cint, Cint, Cint,
               Ptr{Float64}),
              ptr(hh), ptr(ae), loc === nothing ? C_NULL : ptr(loc), data,
              sizeof(data), l, MODES[mode], tests, ee)
    end
    done!(hh, ifail, check)
end

"""
the local potential alone of the input, from the all-electron atom of `ae`
(oncv_local_toml): fields local.*, for `channel(...; loc)` and `assemble`
(`local` is Julia's keyword, so the argument is `loc`)
"""
function local_potential(inp, ae::Handle; check=true, text=true, tables...)
    data = input_text(inp, tables)
    hh = Local(text)
    ifail = ccall(sym(:oncv_local_toml), Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{UInt8}, Csize_t),
                  ptr(hh), ptr(ae), data, sizeof(data))
    done!(hh, ifail, check)
end

"""
the run of the input from the all-electron atom of `ae` with the channels
made before (`channel`, modes :channel or :projectors, in any order) in
place of optimizing them (oncv_assemble_toml): one of each l from 0 to
lmax, the local channel (lloc < 4) `loc`'s unless one of the channels
is of l = lloc. Each must be of this input and atom; the run is then
`oncvpsp(inp; ae)`'s bit for bit, its text output without the channels'
lines
"""
function assemble(inp, ae::Handle, channels; loc=nothing, skip_ld=false, check=true,
                  text=true, tables...)
    data = input_text(inp, tables)
    ptrs = Ptr{Cvoid}[ptr(ch) for ch in channels]
    hh = Run(text)
    ifail = GC.@preserve channels ccall(sym(:oncv_assemble_toml), Cint,
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}, Cint, Ptr{UInt8}, Csize_t,
         Cint),
        ptr(hh), ptr(ae), loc === nothing ? C_NULL : ptr(loc), ptrs, length(ptrs),
        data, sizeof(data), skip_ld)
    done!(hh, ifail, check)
end

"""
whether the first pseudo-wave of channel `l` has a node (pspot's test),
from the all-electron atom of `ae`; a failure before the test throws
"""
function has_node(inp, ae::Handle, l::Integer; text=false, tables...)
    hh = channel(inp, ae, l; mode=:node, check=false, text=text, tables...)
    hh.ifail == 0 && return false
    err = error_of(hh)
    isnode(err) && return true
    throw(err)
end

"""
the model core, screening and tests again on the stages of `base`'s run
with the model core icmod, fcfact, rcfact and r_s (icmod 8, 9; 0 for the
input's) (oncv_coupled): a run
"""
function coupled(base::Run, icmod::Integer; fcfact=0.0, rcfact=0.0, r_s=0.0,
                 skip_ld=false, check=true, text=true)
    hh = Run(text)
    ifail = ccall(sym(:oncv_coupled), Cint,
                  (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Cdouble, Cdouble, Cdouble, Cint),
                  ptr(hh), ptr(base), icmod, fcfact, rcfact, r_s, skip_ld)
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
function confine(base::Run; radii=nothing, channels=nothing, energies=nothing,
                 walls=nothing, r0=0.0, ecutmax=0.0, check=true, text=true)
    rr = radii === nothing ? Float64[] : Float64[radii...]
    nradii = radii === nothing ? 0 : (isempty(rr) ? -1 : length(rr))
    cc = channels === nothing ? Cint[] : Cint[channels...]
    energies === nothing || length(energies) == 3 ||
        throw(ArgumentError("energies is (emin, emax, de)"))
    ee = energies === nothing ? C_NULL : Float64[energies...]
    ww = walls === nothing ? Float64[] : Float64[walls...]
    nwalls = walls === nothing ? -1 : length(ww)
    hh = Confinement(text)
    ifail = GC.@preserve rr cc ee ww ccall(sym(:oncv_confine), Cint,
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Float64}, Cint, Ptr{Cint}, Cint, Ptr{Float64},
         Ptr{Float64}, Cint, Cdouble, Cdouble),
        ptr(hh), ptr(base), rr, nradii, cc, length(cc), ee, ww, nwalls, r0, ecutmax)
    done!(hh, ifail, check)
end

end # module
