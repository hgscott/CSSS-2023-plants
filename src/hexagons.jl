using Hexagons,Colors,Plots

include("Cell.jl")

model = Model(01,Set(),Set())
#dt,cells,nodes

p = plot(legend = false)
for i in 1:40
    for j in 1:4
        hex = HexagonAxial(i, j)  # Creates a hexagonal grid of size 3
        verti = vertices(hex)
        cell_nodes = []
        for vert in verti
            ids = [n.id for n in collect(model.nodes)]
            if length(ids) == 0
                id = -1
            else
                id = minimum(ids) - 1
            end
            bool,node = get_node(vert,model)
            if !bool
                node = Node(id,vert,Set(),Dict())
                push!(model.nodes,node)
            end
            push!(cell_nodes,node)
        end
        ids = [c.id for c in collect(model.cells)]
        if length(ids) == 0
            id = 0
        else
            id = maximum(ids) + 1
        end
        cell = Cell(id,(0,0),1 + 10 * j / 4,0,
        cell_nodes,
        Dict())
        push!(model.cells,cell)
        sync_partners!(cell,model)
    end
end

plot_cells(Plots.plot(legend=false),model)


for i = 1:5_000
    if i % 100 == 0
        p = plot_cells(Plots.plot(legend=false),model)
        display(p)
        display(i)
    end
    relax!(model)
end

plot_cells(Plots.plot(legend=false),model,true)

(c -> c.dV).(collect(model.cells))



loadpath = "./plots/"
# path to where plots should be stored
anim = Animation(loadpath,String[])

for i = 1:20_000
    relax!(model)
    # relax_nodes!(model)
    p = plot_cells(Plots.plot(legend=false),model)
    if i % 80 == 0
        display(i)
        frame(anim,p)
    end
end

gif(anim)

vid_name = "test.mp4"
run(`ffmpeg -framerate 15 -i $loadpath"%06d.png" -vcodec libx264  -crf 25 $vid_name`)
