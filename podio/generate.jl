using YAML
using Graphs

const builtin_types = Dict("int" => "Int32", "float" => "Float32", "double" => "Float64",
    "bool" => "Bool", "long" => "Int64", "unsigned int" => "UInt32", 
    "uint32_t" => "UInt32", "int32_t" => "Int32",  
    "uint64_t" => "UInt64", "int64_t" => "Int64",
    "uint16_t" => "UInt16", "int16_t" => "Int16",
    "uint8_t" => "UInt8", "int8_t" => "Int8",
    "unsigned long" => "UInt64", "char" => "Char", "short" => "Int16",
    "long long" => "Int64", "unsigned long long" => "UInt64")

# reserved words in Julia that cannot be used as field names
const reserved_words = Set([
    "break", "catch", "continue", "do", "else", "elseif", "end", "for", "function",
    "if", "import", "in", "let", "local", "macro", "quote", "return", "try", "using",
    "while", "with", "where"
])

const fundamental_types = [
    "Int8", "UInt8", "Int16", "UInt16", "Int32", "UInt32",
    "Int64", "UInt64", "Int128", "UInt128",
    "Float16", "Float32", "Float64",
    "Complex{T}", "Bool", "Char", "Void",
    "Ptr{T}", "Function", "Function{T}",
    "Csize_t", "Cptrdiff_t"
]

const interface_subtypes = Dict{String,String}()

function to_julia(ctype)
    ctype = ctype |> strip 
    #---Primitive type
    haskey(builtin_types, ctype) && return builtin_types[ctype]
    #---edm4hep type
    m = match(r"^(edm4hep|edm4eic)::(.*)", ctype)
    !isnothing(m) && return m.captures[1] * "!" * m.captures[2]
    #---std::array type
    m = match(r"std::array<([^,]+)[, ]+([0-9]+)>", ctype)
    !isnothing(m) && return "SVector{$(m.captures[2]),$(to_julia(m.captures[1]))}"
    #---Error
    error("Type [$ctype] not translatable to Julia")
end

function gen_member(v,t)
    vt = "$(v)::$(t)"
    vt = vt * " "^(length(vt) > 32 ? 1 : 32 - length(vt))
end

function split_member(member)
    comment = ""
    units = ""
    defvalue = ""
    # extract comment
    pos = findfirst("//", member)
    if !isnothing(pos)
        comment =  "# " * member[pos[1]+2:end]
        member = member[1:pos[1]-1] |> strip
    end
    # extract default value
    m = match(r"(.*){(.+)}", member)
    if !isnothing(m)
        member = m.captures[1] |> strip
        defvalue = m.captures[2] |> strip
    end
    # extract units     
    m = match(r"(.*)[ ]*(\[.*\])", member)
    if !isnothing(m)
        member = m.captures[1] |> strip
        units = m.captures[2]
    end
    sep = findlast(' ', member)
    var = member[sep+1:end]
    if var in reserved_words
        var = var * "_"
        comment = comment * " Renamed from $(member[sep+1:end]) due to reserved word"
    end
    typ = member[1:sep-1] |> to_julia
    defvalue = isempty(defvalue) ? "zero($typ)" : defvalue
    return typ, var, comment, defvalue
end

function gen_interface(io, key, body)
    jtype = to_julia(key)
    gen_docstring(io, key, body)
    println(io, "abstract type $jtype <: POD")
    println(io, "end\n")
end

function gen_component(io, key, body)
    jtype = to_julia(key)
    gen_docstring(io, key, body)
    println(io, "struct $jtype <: POD")
    members = []
    types = []
    defvalues = []
    for m in body["Members"]
        t, v, c ,d = split_member(m)
        vt = "$(v)::$(t)"
        vt = vt * " "^(32 - length(vt))
        println(io, "    $(vt) $(c)")
        push!(members,v)
        push!(types, t)
        push!(defvalues, d)
    end
    args = join(members, ", ")
    defs = join(["$m=$d" for (m,d) in zip(members, defvalues)], ", ")
    println(io, "    $jtype($(defs)) = new($args)")
    println(io, "end\n")
    # add the converters here
    println(io, "Base.convert(::Type{$(jtype)}, t::Tuple) = $(jtype)(t...)")
    ntype = "@NamedTuple{$(join(["$m::$t" for (m,t) in zip(members,types)],", "))}"
    ninit = join(["v.$m" for m in members],", ")
    println(io, "Base.convert(::Type{$(jtype)}, v::$(ntype)) = $(jtype)($(ninit))\n")
end

function gen_datatype(io, key, dtype)
    jtype = to_julia(key)
    gen_docstring(io, key, dtype)
    if haskey(interface_subtypes, jtype)
        println(io, "struct $jtype <: $(interface_subtypes[jtype])")
    else
        println(io, "struct $jtype <: POD")
    end
    vt = gen_member("index", "ObjectID{$jtype}")
    println(io, "    $(vt) # ObjectID of himself")
    println(io, "    #---Data Members")
    members = []
    defvalues = []
    for m in dtype["Members"]
        t, v, c, d = split_member(m)
        println(io, "    $(gen_member(v,t)) $(c)")
        push!(members,v)
        push!(defvalues, t in fundamental_types ? "$d" : contains(t,"SVector") ? "$d" : t*"()")
    end
    vectormembers = @NamedTuple{varname::String, totype::String}[]
    if haskey(dtype, "VectorMembers")
        println(io, "    #---VectorMembers")
        for (i,r) in enumerate(dtype["VectorMembers"])
            t, v, c = split_member(r)
            vt = gen_member(v, "PVector{$(jtype),$(t),$(i)}")
            println(io, "    $(vt) $(c)")
            push!(members, v)
            push!(defvalues, "PVector{$(jtype),$(t),$(i)}()")
            push!(vectormembers, (varname=v,totype=t))
        end
    end
    relations1toN = @NamedTuple{varname::String, totype::String}[]
    if haskey(dtype, "OneToManyRelations")
        println(io, "    #---OneToManyRelations")
        for (i,r) in enumerate(dtype["OneToManyRelations"])
            t, v, c = split_member(r)
            vt = gen_member(v, "Relation{$(jtype),$(t),$(i)}")
            println(io, "    $(vt) $(c)")
            push!(members, v)
            push!(defvalues, "Relation{$(jtype),$(t),$(i)}()")
            push!(relations1toN, (varname=v, totype=t))
        end
    end
    relations1to1 = @NamedTuple{varname::String, totype::String}[]
    if haskey(dtype, "OneToOneRelations")
        println(io, "    #---OneToOneRelations")
        for r in dtype["OneToOneRelations"]
            t, v, c = split_member(r)
            vt = gen_member(v*"_idx", "ObjectID{$(t)}")
            println(io, "    $(vt) $(c)")
            push!(members, v)
            push!(defvalues, "-1")
            push!(relations1to1, (varname=v, totype=t))
        end
    end
    println(io, "end\n")

    # add an extra constructor with keyword parameters
    args = join(members, ", ")
    defs = join(["$m=$dv" for (m,dv) in zip(members,defvalues)], ", ")
    println(io, """
                function $jtype(;$(defs))
                    $jtype(-1, $args)
                end
                """)
    # add an Base.getproperty() for the one-to-one relations
    if !isempty(relations1to1)
        println(io, "function Base.getproperty(obj::$jtype, sym::Symbol)")
        for (i, r) in enumerate(relations1to1)
            if i == 1
                println(io, "    if sym == :$(r.varname)")
            else
                println(io, "    elseif sym == :$(r.varname)")
            end
            println(io, "        idx = getfield(obj, :$(r.varname)_idx)")
            println(io, "        return iszero(idx) ? nothing : convert($(r.totype), idx)")
        end
        println(io, """
                        else # fallback to getfield
                            return getfield(obj, sym)
                        end
                    end""")
    end
    # add pushToXxxx() and popFromXxxx for al one-to-many relations
    if !isempty(relations1toN)
        global exports
        for r in relations1toN
            (;varname, totype) = r
            upvarname = uppercasefirst(varname)
            println(io, "function pushTo$(upvarname)(c::$jtype, o::$totype)")
            println(io, "    iszero(c.index) && (c = register(c))")
            println(io, "    c = @set c.$(varname) = push(c.$varname, o)")
            println(io, "    update(c)")
            println(io, "end")
            println(io, "function popFrom$(upvarname)(c::$jtype)")
            println(io, "    iszero(c.index) && (c = register(c))")
            println(io, "    c = @set c.$(varname) = pop(c.$varname)")
            println(io, "    update(c)")
            println(io, "end")
            push!(exports, "pushTo$(upvarname)", "popFrom$(upvarname)")
        end
    end
    if !isempty(vectormembers)
        global exports
        for v in vectormembers
            (;varname, totype) = v
            upvarname = uppercasefirst(varname)
            println(io, "function set$(upvarname)(o::$jtype, v::AbstractVector{$totype})")
            println(io, "    iszero(o.index) && (o = register(o))")
            println(io, "    o = @set o.$(varname) = v")
            println(io, "    update(o)")
            println(io,"end")
            push!(exports, "set$(upvarname)")
        end
    end
end

function gen_docstring(io, key, dtype)
    jtype = to_julia(key)
    desc = Base.get(dtype, "Description", "$jtype")
    author = Base.get(dtype, "Author", "")
    println(io, "\"\"\"")
    println(io, "$desc")
    !isempty(author) && println(io, "- Author: $author")
    println(io, "# Fields")
    for m in dtype["Members"] 
        t, v, c = split_member(m)
        println(io, "- `$v::$t`: $(c[3:end])")
    end
    for m in Base.get(dtype,"VectorMembers", []) 
        t, v, c = split_member(m)
        t = "PVector{$(t)}"
        println(io, "- `$v::$t`: $(c[3:end])")
    end
    if "OneToOneRelations" in keys(dtype) || "OneToManyRelations" in keys(dtype)
        println(io, "# Relations")
        for m in vcat(Base.get(dtype,"OneToOneRelations",[]),Base.get(dtype,"OneToManyRelations",[])) 
            t, v, c = split_member(m)
            println(io, "- `$v::$t`: $(c[3:end])")
        end
    end
    if !isempty(intersect(("VectorMembers", "OneToManyRelations"),keys(dtype)))
        println(io, "# Methods")
        for m in Base.get(dtype,"VectorMembers", [])
            t, v, c = split_member(m)
            println(io, "- `set$(uppercasefirst(v))(object::$jtype, v::AbstractVector{$t})`: assign a set of values to the `$v` vector member")
        end
        for m in Base.get(dtype,"OneToManyRelations", [])
            t, v, c = split_member(m)
            println(io, "- `pushTo$(uppercasefirst(v))(obj::$jtype, robj::$t)`: push related object to the `$v` relation")
            println(io, "- `popFrom$(uppercasefirst(v))(obj::$jtype)`: pop last related object from `$v` relation")
        end
    end
    println(io,"\"\"\"")
end

function gen_alias(io, dtypes)
    global exports
    jtypes = [to_julia(k) for k in dtypes]
    ns, cs  = zip([split(t,"!") for t in jtypes]...)
    println(io, "\n#---Aliases for easier access")
    for ct in unique(cs)
        ct == "ObjectID" && continue
        ind = findall(x -> x == ct, cs)
        if length(ind) > 1   # In case of name clash, use edm4eic!Xxx
            println(io, "const $(cs[ind[1]]) = edm4eic!$(ct)")
        else
            println(io, "const $(cs[ind[1]]) = $(jtypes[ind[1]])")
        end
        push!(exports, cs[ind[1]])
    end
    println(io, "")
end

function build_graph(datatypes, interfaces=Dict())
    types = to_julia.(keys(datatypes))
    interfaces = to_julia.(keys(interfaces))
    graph = SimpleDiGraph(length(types))
    for (i,dtype) in enumerate(values(datatypes))
        for r in [get(dtype,"OneToOneRelations",[]);get(dtype,"OneToManyRelations",[]);get(dtype,"Members",[])]
            t = split_member(r)[1]
            t == "POD" && continue
            t in interfaces && continue
            d = findfirst(x->x == t, types)
            !isnothing(d) && i != d && add_edge!(graph, d, i)
        end
    end
    # Detect cycles
    scc = strongly_connected_components(graph)
    for (i, dtype) in enumerate(datatypes)
        if any(i in comp for comp in scc if length(comp) > 1)
            println("Datatype [$dtype] is part of a cycle")
        end
    end
    return graph
end

exports = []

function gen_model(sections=[])
    global exports
    #---Load YAML files-------------------------------------------------------------------------------
    yamls = [YAML.load_file(joinpath(@__DIR__, "edm4$(section).yaml")) for section in sections]
    ident = sections[end]

    #---Merge-----------------------------------------------------------------------------------------
    data = Dict()
    data["schema_version"] = yamls[end]["schema_version"]                  # take the schema version of the last one
    data["components"] = merge(getindex.(yamls,"components")...)
    data["datatypes"]  = merge(getindex.(yamls,"datatypes")...)
    data["options"]    = merge(getindex.(yamls,"options")...)
    data["interfaces"] = merge(getindex.(yamls,"interfaces")...)
    data["links"]      = merge(getindex.(yamls,"links")...)

    #---Components-------------------------------------------------------------------------------------
    io = open(joinpath(@__DIR__, "genComponents_$(ident).jl"), "w")

    schema_version = data["schema_version"]
    println(io, "# Automatically generated by generate.jl from edm4hep.yaml (schema version $schema_version)")
    println(io, "schema_version = v\"$schema_version\"\n")

    components = data["components"]
    exports = []
    ctypes = collect(keys(components))
    graph = build_graph(components)
    for i in topological_sort(graph)
        gen_component(io, ctypes[i], components[ctypes[i]])
        push!(exports, to_julia(ctypes[i])) 
    end
    gen_alias(io, ctypes)
    println(io, "\n#---Exports")
    println(io, "export $(join(exports,", "))")
    close(io)

    #---Interfaces-------------------------------------------------------------------------------------
    io = open(joinpath(@__DIR__, "genInterfaces_$(ident).jl"), "w")
    interfaces = data["interfaces"]
    exports = []
    for (key,body) in pairs(interfaces)
        gen_interface(io, key, body)
        push!(exports, to_julia(key))
        for t in body["Types"]
            interface_subtypes[to_julia(t)] = to_julia(key)
        end
    end
    println(io, "\n#---Exports")
    println(io, "export $(join(unique(exports),", "))")
    close(io)

    #---Datatypes--------------------------------------------------------------------------------------
    io = open(joinpath(@__DIR__, "genDatatypes_$(ident).jl"), "w")
    datatypes = data["datatypes"]
    exports = []
    dtypes = collect(keys(datatypes))
    graph = build_graph(datatypes)
    for i in topological_sort(graph)
        gen_datatype(io, dtypes[i], datatypes[dtypes[i]])
    end
    gen_alias(io, dtypes)
    println(io, "\n#---Exports")
    println(io, "export $(join(unique(exports),", "))")
    close(io)

    #---Links------------------------------------------------------------------------------------------
    io = open(joinpath(@__DIR__, "genLinks_$(ident).jl"), "w")
    println(io, "# Automatically generated by generate.jl from edm4hep.yaml (schema version $schema_version)\n")
    links = data["links"]
    exports = []
    for (key,body) in pairs(links)
        jtype = to_julia(key)
        from = to_julia(body["From"])
        to = to_julia(body["To"])
        println(io, "const $jtype = Link{$from,$to}")
    end
    gen_alias(io, keys(links))
    println(io, "\n#---Exports")
    println(io, "export $(join(unique(exports),", "))")
    close(io)
end

gen_model(["hep"])
#gen_model(["hep", "eic"])