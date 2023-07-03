using Agents, LinearAlgebra
using Random,ModularIndices

# both agents need the same number of fields
# so we cann them both at the same time to the ABModel
# maybe create once, then just "subclass" different
# agent types
# dynamics happens in agent_step! anyhow

mutable struct Node <: AbstractAgent
    id::Int
    # :node
    pos::NTuple{2,Float64}
    partners::Set{AbstractAgent} 
    # this should be Cells
end

mutable struct Cell <: AbstractAgent
    id::Int
    # (needed by Agentsjl but not abs. necessary)
    pos::NTuple{2,Float64}

    V::Float64
    dV::Float64


    # name partners more meaningfull
    partners::Vector{Node}
    forces::Dict{AbstractAgent,NTuple{2,Float64}}
    # this should be Nodes
    # both because partners are ordered (in clock direction )

    # if we ever simulate the "walls" as springs. 
    # we need a new class, which "Walls" which goes here
    # again: shared vector of Wall-agents?

end

# maybe add function which adds adds a cell, then adds the nodes
# then pushes the nodes as partners to the cell and vice-versa


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

function calculate_polygon_area(vertices)
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

function agent_step!(cell::Cell, model::ABM)
    # somehow have to calculate dV as 
    # from the node movements
    # probably ned a node.d_pos variable
    

    nodes = cell.partners
    n = length(nodes)
    cell.dV = cell.V - calculate_polygon_area((n -> n.pos).(nodes)) 
    # 10 # * agent.growthrate * model.hardness

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



function split!(i,j,cell,model)
    nodes = cell.partners
    n = length(nodes)
    ii = 1:n 

    function add_node(ij,nodes,model)
        if length(ij) > 1
            pos = .5 .* (nodes[ij[1]].pos .+ nodes[ij[2]].pos)
            n_ids = [n.id for n in allagents(model) if n isa Node]
            id = minimum(n_ids)-1
            n = Node(id,pos,Set())
            add_agent!(n,model)
            if ij[2] < ij[1]
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
    cell_j_partners = ci(nodes,j[end],i[1])

    m_id = maximum([c.id for c in allagents(model) if c isa Cell])

    cell_i = Cell(m_id+1,cell_i_partners[1].pos,2,0,
                cell_i_partners,
                Dict(n => (0,0) for n in cell_i_partners))
    cell_j = Cell(m_id+2,cell_j_partners[1].pos,2,0,
                cell_i_partners,
                Dict(n => (0,0) for n in cell_j_partners))

    
    # still need to change this
    # such that 


    push!(node1.partners,cell_i)
    push!(node2.partners,cell_j)

    
    add_agent_pos!(cell_i,model)
    add_agent_pos!(cell_j,model)
        

    remove_agent!(cell, model)


end


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

    cell = cells[1]
    split!(1,(2,3),cell,model)



    for c in cells
        if c.dV < .001
            split!(1,(2,3),c,model)
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


space = ContinuousSpace((4,4), spacing = 1.0; periodic = false)
model = ABM(
    Union{Cell,Node},
    space,
    properties = Dict(:dt => 0.001, :hardness => 1e2, :mobility => 1.0),
    rng = MersenneTwister(1680)
)



n1 = Node(-1,(1,1),Set())
n2 = Node(-2,(1,2),Set())
n3 = Node(-3,(2,1),Set())


c1 = Cell(1,(0,0),2,0,
            [n1,n2,n3],
            Dict(n1=>(0,0), n2=> (0,0), n3=>(0,0)))


push!(n1.partners,c1)
push!(n2.partners,c1)
push!(n3.partners,c1)

add_agent_pos!( c1,model)
add_agent_pos!( n1,model)
add_agent_pos!( n2,model)
add_agent_pos!( n3,model)


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

function get_poly(agent)
    if agent isa Node
        return Nothing
    elseif agent isa Cell
        coords = [Point2f(p.pos) for p in c1.partners]
        return Polygon(coords)
    end
end

function plot_model(model)
    cells = [c for c in allagents(model) if c isa Cell]
    for cell in cells
        poly!(get_poly(cell), 
                color = :red, strokecolor = :black, strokewidth = 1)
    end
end


fig = Figure()
ax = Axis(fig[1, 1])

nframes = 1000
framerate = 30
it = range(0, nframes)

record(fig, "color_animation.mp4", it;
        framerate = framerate) do _
    plot_model(model)
    model_step!(model)
end
