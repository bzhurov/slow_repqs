NDIM = 19                          #system dimension
IPS = 1                            #this constant defines the problem type
IRS = 0                            #This constant sets the label of the solution where the computation is to be restarted
ILP = 0                            #detect folds
ICP = [ 1, 2 ]                     #free parameters
NTST = 500                         #The number of mesh intervals used for discretization
NCOL = 4                           #The number of Gauss collocation points per mesh interval
IAD = 3                            #Adapt the mesh every IAD
ISP = 0                            #This setting enables the detection of all special solutions
ISW = 2                            #branch switching - for codim 1 bifurcations
IPLT = 0                           #L2 norm will be printed
NBC = 0                            #The number of boundary conditions as specified in the user-supplied routine BCND
NINT = 0                           #The number of integral conditions as specified in the user-supplied routine ICND
LAB = 1
IBR = 1

NMX = 1000000                      #The maximum number of steps to be taken along any family
NPR = 2                            #This constant can be used to regularly write fort.8
MXBF = 10                          #This constant, which is effective for algebraic problems only, sets the maximum number of bifurcations to be treated
IID = 0                            #NO!! diagnostic output to fort.9
ITMX = 100                         #The maximum number of iterations allowed in the accurate location of special solutions
ITNW = 5                           #The maximum number of combined Newton-Chord iterations
NWTN = 3                           #After NWTN Newton iterations the Jacobian is frozen
JAC = 1                            #Derivatives with respect to state- and problem-parameters are given in the usersupplied routines

EPSL = 1e-05                       #Relative convergence criterion for equation parameters in the Newton/Chord method
EPSU = 1e-05                       #Relative convergence criterion for solution components in the Newton/Chord method
EPSS = 1e-06                       #Relative arclength convergence criterion for the detection of special solutions

DS = 0.1                           #
DSMIN = 1e-20                      #This is minimum allowable absolute value of the pseudo-arclength stepsize
DSMAX = 2                          #The maximum allowable absolute value of the pseudo-arclength stepsize
IADS = 1                           #Adapt the pseudo-arclength stepsize after every IADS steps

NPAR= 4                            #
THL =  {}                          #Weights in step-size definition
THU =  {}                          #
UZSTOP = {1: 0, 2: 0}
