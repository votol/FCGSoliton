This project is a calculation programm to numericaly investigate 
qunatum solitons present in onedimensional array of singlemode  
waveguides with cubic nonlinearity. The evolution of quantum properties  
of light is described via a set of nonlinear ordinary differential equations
obtaind in Gaussion approximation for quantum state of light. 

To build code use 
    
    make
    
To run the result
    
    make run_test
    
or better use 
    
    make release
    
after which in workDir will appear a directory with name revN, where N is some number form 1.
In this directory will appear binary called CalcProj, it requaires .yaml file as it only input.
See "test" folder for example, also please see description.yaml file for description of what should be placed 
in "parameters" and "properties" section of this input .yaml file.
