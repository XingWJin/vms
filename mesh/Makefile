all: square geoms meshes pumi_meshes

square:
	mpicxx -std=c++11 square.cpp -o square -I/Users/bng/gmodel/install-opt/include -L/Users/bng/gmodel/install-opt/lib -lgmodel

geoms:
	./square 10
	./square 20
	./square 40
	./square 80
	./square 160

meshes:
	gmsh -2 square_10.geo
	gmsh -2 square_20.geo
	gmsh -2 square_40.geo
	gmsh -2 square_80.geo
	gmsh -2 square_160.geo

pumi_meshes:
	from_gmsh square_10.dmg square_10.msh square_10_.smb
	from_gmsh square_20.dmg square_20.msh square_20_.smb
	from_gmsh square_40.dmg square_40.msh square_40_.smb
	from_gmsh square_80.dmg square_80.msh square_80_.smb
	from_gmsh square_160.dmg square_160.msh square_160_.smb

clean:
	rm -f *.dmg *.geo *.msh square *.smb
