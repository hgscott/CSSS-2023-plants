using LinearAlgebra,Statistics
using Random,ModularIndices,Plots


# both agents need the same number of fields
# so we cann them both at the same time to the ABModel
# maybe create once, then just "subclass" different
# agent types
# dynamics happens in agent_step! anyhow

mutable struct Node
    id::Int
    # :node
    pos::NTuple{2,Float64}
    partners::Set{Int}
    forces::Dict{Int,NTuple{2,Float64}}
    # this should be Cells
end

mutable struct Cell
    id::Int
    # (needed by Agentsjl but not abs. necessary)
    pos::NTuple{2,Float64}

    V::Float64
    dV::Float64


    # name partners more meaningfull
    # here i am using circular refences
    # should probalby use ids instead
    # and then define a "partners" function
    # using multiple dispatch
    partners::Vector{Node}
    forces::Dict{Node,NTuple{2,Float64}}
    # this should be Nodes
    # both because partners are ordered (in clock direction )

    # if we ever simulate the "walls" as springs. 
    # we need a new class, which "Walls" which goes here
    # again: shared vector of Wall-agents?

end

mutable struct Model
    dt::Float64
    cells::Set{Cell}
    nodes::Set{Node}

end

function get_cell(id,model::Model)
    for cell in collect(model.cells)
        if cell.id == id
            return cell
        end
    end
    error("Cell-id not in model")
end

function get_node(id,model::Model)
    for node in collect(model.nodes)
        if node.id == id
            return node
        end
    end
    error("Node-id not in model")
end

function get_partners(node::Node,model::Model)
   return Set( (id -> get_cell(id,model)).(node.partners) )
end

# maybe add function which adds adds a cell, then adds the nodes
# then pushes the nodes as partners to the cell and vice-versa

function get_h(x::Tuple{Float64,Float64})
    normal_xy = (-x[2],x[1])
    return normal_xy ./ norm(x)
end


function agent_step!(node::Node,model::Model)
    # first step using current fs 
    # (since they are already) adjusted by the model
    f = (0,0)
    for id in node.partners
        f = f .+ get_cell(id,model).forces[node]
    end
    
    node.pos = node.pos .+ model.dt .* f

    # here we have to change dV of 
    # the partner cells
    return node.pos
end

function circ_fx(f,x)
    n = length(x)
    res = Vector(undef,n)
    for i in 1:(n-1)
        res[i] = f(x[i],x[i+1])
    end
    res[n] = f(x[n],x[1])
    return res
end

function calculate_polygon_area(vertices::Vector{NTuple{2,Float64}})
    n = length(vertices)
    area = 0.0

    for i in 1:n
        x1, y1 = vertices[i]
        x2, y2 = vertices[mod(i % n, n) + 1]  # Wrap around to the first vertex

        area += (x1 * y2 - x2 * y1)
    end

    area /= 2.0
    return abs(area)
end

function agent_step!(cell::Cell, model::Model)
    # somehow have to calculate dV as 
    # from the node movements
    # probably ned a node.d_pos variable

    nodes = cell.partners
    n = length(nodes)
    cell.dV = cell.V - calculate_polygon_area((n -> n.pos).(nodes)) 
    # 10 # * agent.growthrate * model.hardness

    # connecting vector
    x_i = circ_fx((n1,n2) -> (n2.pos .- n1.pos),nodes)
    # length of sides
    nx_i = norm.(x_i)

    x0_i = (x -> x[1] ./ x[2]).(zip(x_i,nx_i))

    # strength due to "pressure" (dV) per wall
    dh = cell.dV / sum(nx_i)

    h_i = get_h.(x0_i)
    fh_i = (i -> nx_i[i] .* h_i[i]).(1:n)

    forces = circ_fx((fh_1,fh_2) -> (fh_1 .+ fh_2),circshift(fh_i,1))

    forces = (f -> f ./ sum(norm.(forces))).(forces)

    # calculating spring force 
    rel_nx_i =  mean(nx_i)  .- nx_i
    spring_forces = (x -> x[1] .* x[2]).(zip(rel_nx_i,x0_i))
    spring_forces = circ_fx((f_prev,f_next) -> (f_prev .- f_next),circshift(spring_forces,1))
    spring_forces = (f -> f ./ sum(norm.(spring_forces))).(spring_forces)

    forces = (x -> (1+cell.dV) .* x[1] .+ x[2]).(zip(forces,spring_forces))


    for i in 1:n
        cell.forces[nodes[i]] = forces[i]
    end
    cell.pos = cell.partners[1].pos


    return cell.pos
end


function sync_partners!(cell::Cell,model::Model)
    for p in cell.partners
        display(p.id)
        push!(p.partners,cell.id)
    end
end



function split!(i::Tuple{Int,Int},j::Tuple{Int,Int},
                cell::Cell,model::Model)
    nodes = cell.partners
    n = length(nodes)

    i_partnercell = intersect(
                    get_partners(nodes[i[1]],model),
                    get_partners(nodes[i[2]],model)
                                    )
    j_partnercell = intersect(
                    get_partners(nodes[j[1]],model),
                    get_partners(nodes[j[2]],model)
                                    )

    delete!(i_partnercell,cell)
    delete!(j_partnercell,cell)

    ii = 1:n
    i,j = sort([j,i])
    j = j .+ 1

    function add_node(ij,nodes,model)
        if length(ij) > 1
            pos = .5 .* (nodes[ij[1]].pos .+ nodes[ij[2]].pos)
            n_ids = [n.id for n in model.nodes if n isa Node]
            id = minimum(n_ids)-1
            n = Node(id,pos,Set(),Dict())
            push!(model.nodes,n)
            if ij[1] < ij[2]
                return n,vcat(nodes[1:ij[1]],n,nodes[ij[2]:end])
            else
                return n,vcat(nodes,n)
            end
        elseif length(ij) == 1
            n = nodes[ij]
            return n,nodes
        end
    end

    node1,nodes = add_node(i,nodes,model)
    node2,nodes = add_node(j,nodes,model)


    function ci(a,i,j)
        if i < j
            return a[i:j]
        else
            return vcat(a[i:end],a[1:j])
        end
    end

    # if len(i/j) == 1
    # return node
    cell_i_partners = ci(nodes,i[end],j[end])
    cell_j_partners = ci(nodes,j[end],i[end])

    m_id = maximum([c.id for c in model.cells])

    cell_i = Cell(
        m_id+1,cell_i_partners[1].pos,2,0,
                cell_i_partners,Dict())
    cell_j = Cell(
        m_id+2,cell_j_partners[1].pos,2,0,
                cell_j_partners,Dict())
    
    push!(model.cells,cell_i,)
    push!(model.cells,cell_j)

    sync_partners!(cell_i,model)
    sync_partners!(cell_j,model)
        
    
    for n in nodes
        delete!(n.partners,cell.id)
    end
    delete!(model.cells,cell)
end

function update_forces!(model::Model)
    for c in model.cells
        agent_step!(c,model)
    end
end

function relax_nodes!(model::Model)
    #goes through all cell agents and updates forces
    # this can be parallelize
    #this as well
    for n in model.nodes
        # display(n.partners)
        agent_step!(n,model)
    end
end

function relax!(model::Model)
    #goes through all cell agents and updates forces
    # this can be parallelized
    update_forces!(model)
    #this as well
    relax_nodes!(model)
end

function model_step!(model::Model,splitting::Bool)

    relax!(model)

    if splitting
        for c in model.cells
            if c.dV < .001 && Random.rand(1)[1] < .05
                # display("splitting")
                split!((1,2),(2,3),c,model)
            end
        end
    end
    # however, important 
    # all cells then all nodes 
    # or vice versa but NOT: cell1,node1,cell2,node2,node3,....

    # than we need something which goes through all nodes and moves them
    # according to the force
    # Q: HOW do we come up with a Vector of all Nodes.
    # ? Maybe a global Vector (can we change it then?)
    # ...Mutability? has to be updated after a split
end
nothing
####
# creating an example

function cell_mean(cell::Cell)
    nodes = cell.partners
    mean_x = mean([n.pos[1] for n in nodes])
    mean_y = mean([n.pos[2] for n in nodes])
    return mean_x,mean_y
end

# model_step!(model)
function plot_model(model::Model)
    p = Plots.plot()

    nodes = model.nodes
    n = length(nodes)
    nodes_pos = [n.pos for n in nodes]
    n_x = (n->n[1]).(nodes_pos)
    n_y = (n->n[2]).(nodes_pos)
    p = Plots.scatter!(p,n_x,n_y,
            markersize = 0)

    for n in model.nodes
        Plots.annotate!(p,n.pos[1],n.pos[2],
                Plots.text(string(n.id), :black, 20))
    end

    c_color = :green
    for c in model.cells
        Random.seed!(c.id)
        color = rand(RGB)
        m = cell_mean(c)
        Plots.annotate!(p,m[1],m[2],
            Plots.text(string(c.id),:black, 20))
        
        n = length(c.partners)
        nodes_pos = [n.pos for n in c.partners]
        n_x = (n->n[1]).(nodes_pos)
        n_y = (n->n[2]).(nodes_pos)
        Plots.plot!(p,
            n_x[Mod(1:(n+1))],n_y[Mod(1:(n+1))],
            fill = (0, color),)
    end

    return p
end



n1 = Node(-1,(1,1),Set(),Dict())
n2 = Node(-2,(1,10),Set(),Dict())
n3 = Node(-3,(1.1,1),Set(),Dict())


c1 = Cell(1,(0,0),2,0,
            [n1,n2,n3],
            Dict(n1=>(0,0), n2=> (0,0), n3=>(0,0)))


push!(n1.partners,c1.id)
push!(n2.partners,c1.id)
push!(n3.partners,c1.id)

model = Model(0.01,Set([c1]),Set([n1,n2,n3])) #dt,cells,nodes

p = plot_model(model)
display(p)

for i = 1:1000
    relax!(model)
   # relax_nodes!(model)
   if i % 20 == 0
        p = plot_model(model)
        xs = [n.pos[1] for n in collect(model.nodes)]
        ys = [n.pos[2] for n in collect(model.nodes)]

        xlims!( minimum(xs)-1,maximum(xs)+1 )
        ylims!( minimum(ys)-1,maximum(ys)+1 )
        display(p)
   end
    sleep(.01)
end

cells = collect(model.cells)
nodes = collect(model.nodes)
c = shuffle(cells)[1]
split!((1,2),(2,3),c,model)
p = plot_model(model)

model.nodes
# split!((1,2),(2,3),c1,model) #(i,j,cell,model)


using InteractiveDynamics

function cell_polygon(agent)
    if agent isa Node
        return '.'
    elseif agent isa Cell
        abs_coords = [p.pos for p in c1.partners]
        coords = [ c .- agent.pos for c in abs_coords]
        coords = [Point2f(c) for c in coords]
        return scale(Polygon(coords),8.2)
    end
end
nothing # hide

using GLMakie
using Makie.Colors
# animation settings

function get_poly(cell)
    coords = [Point2f(p.pos) for p in cell.partners]
    return Polygon(coords)
end

function plot_model(model)
    for cell in model.cells
        poly!(get_poly(cell), strokecolor = :black, strokewidth = 5,
                transparency = true)
    end
end


fig = Figure()
ax = Axis(fig[1, 1])

nframes = 100
framerate = 5
it = range(0, nframes)

record(fig, "color_animation.mp4", it;
        framerate = framerate) do _

    plot_model(model)
    model_step!(model)
end

p = Plots.plot()
for cell in model.cells
    display(length(cell.partners))
    Plots.scatter!(p,[p.pos .+ .1 .* rand(2) for p in cell.partners],alpha = .2)
end


using GeometryBasics
# Initialize Figure and Axis
fig = Figure(resolution = (800, 600))
ax = Axis(fig[1, 1])

# Initialize polygons
polygon = [Polygon(rand(Point2f, 3)) for _ in 1:10]
colors = [rand(RGB) for _ in 1:10]

# Plot initial polygons
plots = [lines!(ax, p, color = c) for (p, c) in zip(polygon, colors)]

# Start recording
record(fig, "example.mp4", 1:100) do i
    # Update polygons
    for j in 1:10
        polygon[j] = Polygon(rand(Point2f0, 3))
        plots[j][1] = polygon[j]
    end
end