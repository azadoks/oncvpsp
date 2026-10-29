# The Julia package ONCVPSP against the executable, as scripts/python_wrapper.py
# the Python one: the run's text output and psp8 file the executable's byte for
# byte (its output but for the banner); the input as the TOML file, its
# dictionary and tables given as keyword arguments, the same; the arrays against
# the dimensions; the stages against the run; failures thrown or left in the
# handle, and a handle closed.
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
    run = ONCVPSP.run(toml; irel=irel)
    text = run.text
    ONCVPSP.write(run, "lib.psp8")
    psp8 = read("lib.psp8", String)
    at = findfirst(text, exe)
    @test at !== nothing && occursin("build: ", exe[1:first(at)-1])
    # between them, the fully relativistic executable's lines of the
    # spin-orbit projectors' conversion, which the file's writing makes
    between, after = split(exe[last(at)+1:end], "Begin PSPCODE8\n"; limit=2)
    @test startswith(after, psp8)
    @test all(ll -> isempty(ll) || startswith(ll, "l1-1,npr_sr,npr_so"), split(between, '\n'))

    dims = run.dims
    @test dims.version == ONCVPSP.API_VERSION && dims.irel == irel
    @test size(run["pot.rr"]) == (dims.mmax,)
    vkb = run["pot.vkb"]
    @test size(vkb, 1) == dims.mmax && size(vkb, 4) == dims.nkap
    @test run["pot.zion"] isa Float64
    @test !haskey(run, "channel.pswf")
    fields = Dict(nn => run[nn] for nn in keys(run))
    @test length(fields) > 50
    same(hh) = all(hh[nn] == vv for (nn, vv) in fields)

    data = ONCVPSP.TOML.parsefile(toml)
    @test same(ONCVPSP.run(data; irel=irel))
    @test same(ONCVPSP.run(Dict(); irel=irel, (Symbol(kk) => vv for (kk, vv) in data)...))

    ae = all_electron(data; irel=irel)
    @test all(startswith("ae."), keys(ae))
    @test same(ONCVPSP.run(data; ae=ae))
    for ll in 0:dims.lmax
        ch = channel(data, ae, ll; mode=:projectors)
        ll == dims.lloc || @test ch["channel.evkb"][:, 1] == fields["pot.evkb"][:, ll+1, 1]
        @test !has_node(data, ae, ll)
    end
    qcut = data["pseudopotentials"]["qcut"]
    for factor in (2.0, 3.0)
        pp = merge(data["pseudopotentials"], Dict("qcut" => [qcut[1] * factor; qcut[2:end]]))
        ch = channel(data, ae, 0; check=false, pseudopotentials=pp)
        @test has_node(data, ae, 0; pseudopotentials=pp) ==
              (ch.ifail != 0 && ch.failure.msg == ONCVPSP.MSG_PSPOT_NODE)
    end
    @test confine(run; radii=(), walls=())["confine.nmiss"] == -1
    conf = confine(run; walls=())
    @test conf["confine.nmiss"] >= 0 && occursin("\nws_r ", conf.text)

    err = try
        ONCVPSP.run(Dict("oncvpsp" => Dict("atsym" => "O")); irel=irel)
    catch ee
        ee
    end
    @test err isa OncvpspError && ONCVPSP.kind_name(err) == "input" == ONCVPSP.stage_name(err)
    broken = ONCVPSP.run(Dict("oncvpsp" => Dict("atsym" => "O")); irel=irel, check=false)
    @test broken.ifail != 0 && isempty(keys(broken))
    err = try
        confine(run; radii=0.1)
    catch ee
        ee
    end
    @test err isa OncvpspError && ONCVPSP.stage_name(err) == "confine"
    @test occursin("confine_atom: ERROR", sprint(showerror, err))

    close(run)
    @test_throws ArgumentError run["pot.rr"]
    for _ in 1:3
        all_electron(data; irel=irel, text=false)
    end
    GC.gc()
end
