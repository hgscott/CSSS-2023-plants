using Agents, LinearAlgebra
using Random


# both agents need the same number of fields
# so we cann them both at the same time to the ABModel
# maybe create once, then just "subclass" different
# agent types
# dynamics happens in agent_step! anyhow

mutable struct Node <: AbstractAgent
    id::Int
    # :node
    pos::NTuple{2,Float64}
    partners::Vector{AbstractAgent} 
    # this should be Cells
end

mutable struct Cell <: AbstractAgent
    id::Int
    # (needed by Agentsjl but not abs. necessary)
    pos::NTuple{2,Float64}

    dV::Float64


    # name partners more meaningfull
    partners::Vector{AbstractAgent}
    forces::Dict{AbstractAgent,NTuple{2,Float64}}
    # this should be Nodes
    # both because partners are ordered (in clock direction )

    # if we ever simulate the "walls" as springs. 
    # we need a new class, which "Walls" which goes here
    # again: shared vector of Wall-agents?

end

# maybe add function which adds adds a cell, then adds the nodes
# then pushes the nodes as partners to the cell and vice-versa


function model_step!(model)
    #goes through all cell agents and updates forces
    cells = [c for c in allagents(model) if c isa Cell]
    nodes = [n for n in allagents(model) if n isa Node]

    # this can be parallelized
    for c in cells
        agent_step!(c,model)
    end

    #this as well
    for n in nodes
        agent_step!(n,model)
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

function get_h0(x,nx)
    normal_xy = (-x[2],x[1])
    return normal_xy ./ nx
end


function agent_step!(node::Node,model::ABM)
    # first step using current fs 
    # (since they are already) adjusted by the model
    f = (0,0)
    for cell in node.partners
        f = f .+ cell.forces[node]
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

function agent_step!(cell::Cell, model::ABM)
    # somehow have to calculate dV as 
    # from the node movements
    # probably ned a node.d_pos variable
    
    cell.dV = 3 # 10 # * agent.growthrate * model.hardness

    nodes = cell.partners
    n = length(nodes)

    x_i = circ_fx((n1,n2) -> (n2.pos .- n1.pos),nodes)

    nx_i = norm.(x_i)
    dh = cell.dV / sum(nx_i)

    h_i = (tup -> get_h0(tup[1],tup[2])).(zip(x_i,nx_i))
    fh_i = (i -> nx_i[i] .* h_i[i]).(1:n)
    fh_i = circshift(fh_i,1)

    forces = circ_fx((fh_1,fh_2) -> dh .* (fh_1 .+ fh_2),fh_i)

    for i in 1:n
        cell.forces[nodes[i]] = forces[i]
    end

    cell.pos = cell.partners[1].pos


    return cell.pos
end

####
# creating an example


space = ContinuousSpace((4,4), spacing = 1.0; periodic = false)
model = ABM(
    Union{Cell,Node},
    space,
    properties = Dict(:dt => 0.001, :hardness => 1e2, :mobility => 1.0),
    rng = MersenneTwister(1680)
)



n1 = Node(-1,(1,1),[])
n2 = Node(-2,(1,2),[])
n3 = Node(-3,(2,1),[])

c1 = Cell(1,(0,0),0,
            [n1,n2,n3],
            Dict(n1=>(0,0), n2=> (0,0), n3=>(0,0)))


push!(n1.partners,c1)
push!(n2.partners,c1)
push!(n3.partners,c1)

add_agent_pos!( c1,model)
add_agent_pos!( n1,model)
add_agent_pos!( n2,model)
add_agent_pos!( n3,model)


[p.pos for p in c1.partners]
model_step!(model)
[p.pos for p in c1.partners]





using InteractiveDynamics
using CairoMakie # choose plotting backend
CairoMakie.activate!() # hide

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

# set up some nice colors
cell_color(b) = RGBf(b.id * 3.14 % 1, 0.2, 0.2)
nothing # hide

# and proceed with the animation
InteractiveDynamics.abmvideo(
    "cell.mp4", model, agent_step!;
    am = cell_polygon, ac = cell_color,
    spf = 10, framerate = 30, frames = 600,
    title = "Growing plants"
)

