```@meta
EditURL = "tutorial_io.lit"
```

# EDM4hep.jl EDM I/O Tutorial
This tutorial shows how to read EDM4hep data files and perform some simple analysis

!!! note "Note that"
    You can also download this tutorial as a
    [Jupyter notebook](tutorial_io.ipynb) and a plain
    [Julia source file](tutorial_io.jl).

#### Table of contents
```@contents
Pages = ["tutorial_io.md"]
Depth = 2:3
```


#### Loading the necessary modules
- We will use the `EDM4hep` module that defines all the EDM4hep types and its relations and links.
- We will use the [`FHist`](https://juliaphysics.github.io/FHist.jl/stable/) module to create histograms
- We will use the [`Plots`](https://docs.juliaplots.org/stable/) module to plot the histograms
If these modules are not installed, you can install them by running the following commands:
```julia
using Pkg
Pkg.add("EDM4hep")
Pkg.add("FHist")
Pkg.add("Plots")
```

````julia
using EDM4hep
using EDM4hep.RootIO
using Base.Iterators: partition, take
using FHist
using Plots: plot, scatter, plot!, theme
````

## Reading an EDM4hep file
We will read an EDM4hep file with the `RootIO.Reader` function. This function returns a reader object that can be used to access
the events in the file. The input file is a ROOT file with the EDM4hep data model and is located to the path `ttbar_edm4hep_digi.root`.

````julia
finput = joinpath(@__DIR__,"../../examples" ,"ttbar_edm4hep_digi.root")
````

````
"/Users/mato/Development/EDM4hep.jl/docs/src/tutorials/../../examples/ttbar_edm4hep_digi.root"
````

We create a reader object to access the events in the file. The object displays a summary table of the content of the file.

````julia
reader = RootIO.Reader(finput)
````

````


ERROR: SystemError: opening file /Users/mato/Development/EDM4hep.jl/docs/src/tutorials/../../examples/ttbar_edm4hep_digi.root: No such file or directory
````

### Accessing the events
The `TTree` called `events` contains the events in the file. We can access the events by using the `RootIO.get` function
with the reader object and the name of the `TTree`.

````julia
events = RootIO.get(reader, "events");
````

````


ERROR: UndefVarError: `reader` not defined in `Main.var"##277"`
Suggestion: add an appropriate import or assignment. This global was declared but not assigned.
````

!!! note "Note that"
    If you get warnings is because there is a mismatch between the
    current schema version of EDM4hep and version when the file was written.
    The default behavior is ignore old types and set to zero the attributes
    that do not exists in the file.

We can access the first event in the file by using the index `1` in the `events` array.

````julia
evt = events[1];
````

````


ERROR: UndefVarError: `events` not defined in `Main.var"##277"`
Suggestion: add an appropriate import or assignment. This global was declared but not assigned.
````

the `evt` object is a `UnROOT.LazyEvent`. Each leaf (column) can be accessed directlyt if you know the name
of the leaf. Typically the name is of the form `<collection_name>_<field_name>`. The full list of names can be
obtained by calling the `names` function on the `events` object. For example:

````julia
for n in names(events) |> sort!
    startswith(n, "ECalBarrelCollection") && println(n)
end
````

````


ERROR: UndefVarError: `events` not defined in `Main.var"##277"`
Suggestion: add an appropriate import or assignment. This global was declared but not assigned.
````

!!! note "Note that"
    There is not need to access individual columns like this. The `RootIO.get` function can be used to access the collections directly.

### Accessing the collections
The available collections in the event can be obtained by displaying the reader object.

````julia
show(reader)
````

````


ERROR: UndefVarError: `reader` not defined in `Main.var"##277"`
Suggestion: add an appropriate import or assignment. This global was declared but not assigned.
````

The `RootIO.get` function can be used to access the collections in the event.
The function takes the reader object, the event object and the name of the collection as arguments.

````julia
calo = RootIO.get(reader, evt, "AllCaloHitContributionsCombined");
hits = RootIO.get(reader, evt, "ECalBarrelCollection");
mcps = RootIO.get(reader, evt, "MCParticle");
````

````


ERROR: UndefVarError: `reader` not defined in `Main.var"##277"`
Suggestion: add an appropriate import or assignment. This global was declared but not assigned.
````

!!! note "Note that"
    Relationships will automatically be setup if the related collections will be also be accessed.
    For example, in this case the relation of `SimCalorimeterHit`s from the collection `ECalBarrelCollection`
    to `CaloHitContribution`s from `AllCaloHitContributionsCombined` will be filled.

The `calo` object is a `EDM4hep.EDCollection` object that contains the calorimeter hit contributions.
The `hits` object is a `EDM4hep.EDCollection` object that contains the hits in the ECal barrel.
The `mcps` object is a `EDM4hep.EDCollection` object that contains the MCParticles.

We can now print some information about the collections.

````julia
for hit in take(hits, 20)
    println("ECAL Hit $(hit.index) has energy $(hit.energy) at position $(hit.position) with $(length(hit.contributions)) contributions")
end
````

````


ERROR: UndefVarError: `hits` not defined in `Main.var"##277"`
Suggestion: check for spelling errors or missing imports.
````

Lets check whether the total energy in the hits is the same as the sum of all the calorimeter contributions.
This an example to show the expressivity of the Julia language.
We can construct a `StructArray` with all the related hit contributions and use the `sum` function
to sum the energy of all the contributions as a column of the constructed SoA.

````julia
for hit in hits
    StructArray(hit.contributions).energy |> sum |> c -> c ≈ hit.energy || println("Hit $(hit.index) has energy $(hit.energy) and contributions $(c)")
end
````

````


ERROR: UndefVarError: `hits` not defined in `Main.var"##277"`
Suggestion: check for spelling errors or missing imports.
````

### Drawing the EDMM4hep data
The following shows how easy is to draw EDM4hep event using the `Plots` module. In this example we want to plot the calorimeter hits
in the space with dots with the size proportional to the energy of the hit.

````julia
maxE = hits.energy |> maximum
scatter(hits.position.x, hits.position.y, hits.position.z,
        markersize = (hits.energy/maxE)*10, color = :blue)
````

````


ERROR: UndefVarError: `hits` not defined in `Main.var"##277"`
Suggestion: check for spelling errors or missing imports.
````

The projection in the X-Y plane is shown below.

````julia
scatter(hits.position.x, hits.position.y, markersize = (hits.energy/maxE)*10, color = :blue)
````

````


ERROR: UndefVarError: `hits` not defined in `Main.var"##277"`
Suggestion: check for spelling errors or missing imports.
````

Accessing the `hit` attributes as columns of the `hits` object is very efficient.
Lets verify it. Using the `@btime` macro from the `BenchmarkTools` module we can measure the time to access the `energy` column.

````julia
using BenchmarkTools
@benchmark hits.energy |> maximum
````

````


ERROR: ArgumentError: Package BenchmarkTools not found in current path.
- Run `import Pkg; Pkg.add("BenchmarkTools")` to install the BenchmarkTools package.
````

The time to access the `energy` column is very small.

````julia
@benchmark begin
    _max = 0.0f0
    for h in hits
        _max = max(_max, h.energy)
    end
    _max
end
````

````


ERROR: LoadError: UndefVarError: `@benchmark` not defined in `Main.var"##277"`
Suggestion: check for spelling errors or missing imports.
in expression starting at /Users/mato/Development/EDM4hep.jl/docs/src/tutorials/tutorial_io.md:1
````

The time to access the `energy` column using a loop is much larger.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

