all: lshape.smb
lshape.geo:
	mpicxx -std=c++11 lshape.cpp -o lshape -I${GMODEL_DIR}/include -L${GMODEL_DIR}/lib -lgmodel
	./lshape
lshape.msh: lshape.geo
	gmsh -2 lshape.geo -o lshape.msh	
lshape-serial.smb: lshape.msh
	from_gmsh lshape.dmg lshape.msh lshape-serial.smb
lshape.smb: lshape-serial.smb
	mpirun -n 4 split lshape.dmg lshape-serial.smb lshape.smb 4
clean:
	rm -rf lshape lshape.geo lshape.dmg lshape.msh lshape*.smb out*
