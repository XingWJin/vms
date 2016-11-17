NAME1=lshape
NAME2=lshape_square

all: $(NAME1).smb $(NAME2).smb

$(NAME1).geo:
	mpicxx -std=c++11 $(NAME1).cpp -o $(NAME1) -I${GMODEL_DIR}/include -L${GMODEL_DIR}/lib -lgmodel
	./$(NAME1)
$(NAME1).msh: $(NAME1).geo
	gmsh -2 $(NAME1).geo -o $(NAME1).msh  
$(NAME1)-serial.smb: $(NAME1).msh
	from_gmsh $(NAME1).dmg $(NAME1).msh $(NAME1)-serial.smb
$(NAME1).smb: $(NAME1)-serial.smb
	mpirun -n 4 split $(NAME1).dmg $(NAME1)-serial.smb $(NAME1).smb 4

$(NAME2).geo:
	mpicxx -std=c++11 $(NAME2).cpp -o $(NAME2) -I${GMODEL_DIR}/include -L${GMODEL_DIR}/lib -lgmodel
	./$(NAME2)
$(NAME2).msh: $(NAME2).geo
	gmsh -2 $(NAME2).geo -o $(NAME2).msh  
$(NAME2)-serial.smb: $(NAME2).msh
	from_gmsh $(NAME2).dmg $(NAME2).msh $(NAME2)-serial.smb
$(NAME2).smb: $(NAME2)-serial.smb
	mpirun -n 4 split $(NAME2).dmg $(NAME2)-serial.smb $(NAME2).smb 4

clean:
	rm -rf $(NAME1) $(NAME1).geo $(NAME1).dmg $(NAME1).msh $(NAME1)*.smb
	rm -rf $(NAME2) $(NAME2).geo $(NAME2).dmg $(NAME2).msh $(NAME2)*.smb
