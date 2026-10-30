# The Julia package ONCVPSP against the executable, as scripts/python_wrapper.py
# the Python one: the run's text output and psp8 file the executable's byte for
# byte (its output but for the banner); the input as the TOML file, its
# dictionary and tables given as keyword arguments, the same; the arrays against
# the dimensions and each field's axes; the accessors against the fields; the
# stages against the run; failures thrown or left in the handle, and a handle
# closed.
#
# usage: julia --project=julia/ONCVPSP runtests.jl FULL_EXE IREL INPUT.toml WORKDIR
# (ONCVPSP_LIBRARY the library)

using ONCVPSP
using Test

full, irel, toml, work = ARGS[1], parse(Int, ARGS[2]), abspath(ARGS[3]), ARGS[4]
mkpath(work)
cd(work)

Base.write("exe.out", read(`$full --toml-input $toml`))
exe = read("exe.out", String)

@testset "ONCVPSP $(basename(toml)) irel $irel" begin
    run = oncvpsp(toml; irel=irel)
    @test run isa Run && ONCVPSP.run === oncvpsp
    text = run.text
    write("lib.psp8", run)
    psp8 = read("lib.psp8", String)
    at = findfirst(text, exe)
    @test at !== nothing && occursin("build: ", exe[1:first(at)-1])
    # between them, the fully relativistic executable's lines of the
    # spin-orbit projectors' conversion, which the file's writing makes
    between, after = split(exe[last(at)+1:end], "Begin PSPCODE8\n"; limit=2)
    @test startswith(after, psp8)
    @test all(ll -> isempty(ll) || startswith(ll, "l1-1,npr_sr,npr_so"), split(between, '\n'))

    # the library's build, the executable's (one build)
    @test all(occursin("$kk $vv", replace(exe, "," => "")) for (kk, vv) in build_info())

    dims = run.dims
    @test dims.version == ONCVPSP.API_VERSION && dims.irel == irel
    @test size(run["pot.rr"]) == (dims.mmax,)
    vkb = run["pot.vkb"]
    @test size(vkb, 1) == dims.mmax && size(vkb, 4) == dims.nkap
    @test run["pot.zion"] isa Float64
    @test !haskey(run, "channel.pswf")
    fields = Dict(run)
    @test length(fields) > 50 && length(run) == length(fields)
    same(hh) = all(hh[nn] == vv for (nn, vv) in fields)
    axes_ranks(hh) = all(length(fieldinfo(nn).axes) == ndims(vv) && !isempty(fieldinfo(nn).doc)
                         for (nn, vv) in hh)
    @test axes_ranks(run)

    # the accessors against the fields
    @test mesh(run) == fields["pot.rr"]
    @test local_potential(run) == fields["pot.vpuns"][:, dims.lloc+1]
    for ll in 0:dims.lmax, ikap in 1:dims.nkap
        ee, beta = projectors(run, ll; ikap=ikap)
        @test length(ee) == dims.nproj[ll+1]
        kk = sum((ee[ii] * beta[:, ii] * beta[:, ii]' for ii in eachindex(ee));
                 init=zeros(dims.mmax, dims.mmax))
        @test kernel(run, ll; ikap=ikap) ≈ kk
        if !isempty(ee)
            table = convergence(run, ll; ikap=ikap)
            @test size(table, 2) == 2 && all(diff(table[:, 2]) .< 0)
        end
    end
    rr = mesh(run)
    rho = 4π .* rr .^ 2 .* valence_density(run)
    @test abs(sum(diff(rr) .* (rho[1:end-1] .+ rho[2:end]) ./ 2) - fields["pot.zion"]) < 1e-3

    data = ONCVPSP.TOML.parsefile(toml)
    @test same(ONCVPSP.run(data; irel=irel))
    nones = Dict(kk => (vv isa AbstractDict ? merge(vv, Dict("nothing" => nothing)) : vv)
                 for (kk, vv) in data)
    @test same(ONCVPSP.run(nones; irel=irel))
    @test same(ONCVPSP.run(Dict(); irel=irel, (Symbol(kk) => vv for (kk, vv) in data)...))

    ae = all_electron(data; irel=irel)
    @test ae isa AllElectron && all(startswith("ae."), keys(ae)) && axes_ranks(ae)
    @test same(ONCVPSP.run(data; ae=ae))
    for ll in 0:dims.lmax
        ch = channel(data, ae, ll; mode=:projectors)
        @test ch isa ONCVPSP.Channel && axes_ranks(ch)
        if ll != dims.lloc
            @test ch["channel.evkb"][:, 1] == fields["pot.evkb"][:, ll+1, 1]
            @test kernel(ch) ≈ kernel(run, ll)
        end
        @test !has_node(data, ae, ll)
    end
    qcut = data["pseudopotentials"]["qcut"]
    for factor in (2.0, 3.0)
        pp = merge(data["pseudopotentials"], Dict("qcut" => [qcut[1] * factor; qcut[2:end]]))
        ch = channel(data, ae, 0; check=false, pseudopotentials=pp)
        @test has_node(data, ae, 0; pseudopotentials=pp) ==
              (ch.ifail != 0 && ch.failure.msg == ONCVPSP.MSG_PSPOT_NODE)
    end
    # the run assembled from channels made on the local potential alone
    loc = local_potential(data, ae; text=false)
    @test loc isa Local && axes_ranks(loc)
    tiles = ONCVPSP.Channel[]
    for ll in 0:dims.lmax
        ch = channel(data, ae, ll; mode=:projectors, loc=loc, tests=true, text=false)
        mk = size(ch["channel.logder.ae"], 2)
        @test ch["channel.logder.ae"] == fields["tests.logder.ae"][:, ll+1, 1:mk]
        ONCVPSP.isexact(ch) &&
            @test ch["channel.logder.ps"] == fields["tests.logder.ps"][:, ll+1, 1:mk]
        ll == dims.lloc || push!(tiles, ch)
    end
    @test same(assemble(data, ae, tiles; loc=loc, text=false))
    err = try
        assemble(data, ae, tiles[2:end]; text=false)
    catch ee
        ee
    end
    @test err isa OncvpspError && err.kind === :input && err.msg == 190
    @test coupled(run, 0; text=false)["pot.vkb"] == vkb
    @test confine(run; radii=(), walls=())["confine.nmiss"] == -1
    conf = confine(run; walls=())
    @test conf["confine.nmiss"] >= 0 && occursin("\nws_r ", conf.text) && axes_ranks(conf)
    @test oncvpsp(data; irel=irel, text=false) do hh
        hh["pot.vkb"] == vkb
    end

    err = try
        ONCVPSP.run(Dict("oncvpsp" => Dict("atsym" => "O")); irel=irel)
    catch ee
        ee
    end
    @test err isa OncvpspError && err.kind === :input === err.stage && err.msg > 0
    broken = ONCVPSP.run(Dict("oncvpsp" => Dict("atsym" => "O")); irel=irel, check=false)
    @test broken.ifail != 0 && isempty(keys(broken))
    err = try
        confine(run; radii=0.1)
    catch ee
        ee
    end
    @test err isa OncvpspError && err.stage === :confine
    @test occursin("confine_atom: ERROR", sprint(showerror, err))

    close(run)
    @test_throws ArgumentError run["pot.rr"]
    for _ in 1:3
        all_electron(data; irel=irel, text=false)
    end
    GC.gc()
end
