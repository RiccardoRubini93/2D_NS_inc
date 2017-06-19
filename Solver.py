import numpy as np
import math

def build_up_b(b, rho, dt, u, v, dx, dy):
    
    b[1:-1, 1:-1] = (rho*(1/dt * ((u[1:-1, 2:] - u[1:-1, 0:-2])/(2*dx) + (v[2:, 1:-1] - v[0:-2, 1:-1])/(2*dy))-((u[1:-1, 2:] - u[1:-1, 0:-2])/(2*dx))**2-2*((u[2:, 1:-1] - u[0:-2, 1:-1])/(2*dy)*(v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx))-((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2))

    return b


def pressure_poisson(p, dx, dy, b,nit):
    pn = np.empty_like(p)
    pn = p.copy()
    
    for q in range(nit):
        
        pn = p.copy()
        
	#pressure correction technique	
	
        p[1:-1, 1:-1] = (((pn[1:-1, 2:]+pn[1:-1, 0:-2])*dy**2 + (pn[2:, 1:-1]+pn[0:-2, 1:-1])*dx**2)/(2*(dx**2+dy**2))-dx**2 * dy**2/(2*(dx**2+dy**2))* b[1:-1,1:-1])

	# artificial compressibility technique
	#tau = 2 #this parameters has to set un carefully 
	#p[1:-1,1:-1] = pn[1:-1,1:-1] - dt*tau*((u[1:-1,2:] - u[1:-1,0:-2])/(2*dx) + (v[2:, 1:-1] - v[0:-2, 1:-1])/(2*dy))
	
        p[:, -1] = p[:, -2]    ##dp/dy = 0 at x = 2
        p[0, :] = p[1, :]      ##dp/dy = 0 at y = 0
        p[:, 0] = p[:, 1]      ##dp/dx = 0 at x = 0
        p[-1, :] = 0 
    
    return p

def solver(nt, u, v, T, dt, dx, dy, p, rho, nu, alpha,V_top,T_top,nx,ny,nit):
    un = np.empty_like(u)
    vn = np.empty_like(v)
    Tn = np.empty_like(T)
    b = np.zeros((ny, nx))
    
    for n in range(nt):	
	
        resU = np.linalg.norm((un-u)+(vn-v))/math.sqrt(nx*ny)	
        resT = np.linalg.norm(Tn-T)/math.sqrt(nx*ny)	
	
        un = u.copy()
        vn = v.copy()
        Tn = T.copy()
        
        b = build_up_b(b, rho, dt, u, v, dx, dy)
        p = pressure_poisson(p, dx, dy, b,nit)
        
	#Momentum equation solver
	
        u[1:-1,1:-1] = un[1:-1,1:-1] - dt/(2*dx)*(un[1:-1, 2:]*un[1:-1, 2:]-un[1:-1, 0:-2]*un[1:-1, 0:-2]) - dt/(2*dy)*(un[2:, 1:-1]*vn[2:, 1:-1]-un[0:-2, 1:-1]*vn[0:-2, 1:-1]) - dt/(2*rho*dx)*(p[1:-1, 2:]-p[1:-1, 0:-2])+nu*(dt/dx**2*(un[1:-1, 2:]-2*un[1:-1, 1:-1]+un[1:-1, 0:-2])+dt/dy**2*(un[2:, 1:-1]-2*un[1:-1, 1:-1]+un[0:-2, 1:-1]))

        v[1:-1,1:-1] = vn[1:-1,1:-1] -dt/(2*dy)*(vn[1:-1, 2:]*vn[1:-1, 2:]-vn[1:-1, 0:-2]*vn[1:-1, 0:-2]) - dt/(2*dx)*(un[2:, 1:-1]*vn[2:, 1:-1]-un[0:-2, 1:-1]*vn[0:-2, 1:-1]) - dt/(2*rho*dy)*(p[2:, 1:-1]-p[0:-2, 1:-1])+nu * (dt/dx**2 *(vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2])+dt/dy**2 *(vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1]))

	#Energy equation Solver

        T[1:-1,1:-1] = Tn[1:-1,1:-1] -dt/(2*dx)*u[1:-1,1:-1]*(Tn[1:-1, 2:]-Tn[1:-1, 0:-2]) -dt/(2*dy)*v[1:-1,1:-1]*(Tn[2:, 1:-1]-Tn[0:-2, 1:-1]) + alpha*(dt/dx**2*(Tn[1:-1, 2:]-2*Tn[1:-1, 1:-1]+Tn[1:-1, 0:-2]) +dt/dy**2*(Tn[2:, 1:-1] - 2 * Tn[1:-1, 1:-1] + Tn[0:-2, 1:-1]))
	


	#imposing boundary conditions	

        u[0, :] = 0
        u[:, 0] = 0
        u[:, -1] = 0
        u[-1, :] = V_top   
        v[0, :] = 0
        v[-1, :]=0
        v[:, 0] = 0
        v[:, -1] = 0

        T[0, :] = 0
        T[:, 0] = 0
        T[:,-1] = 0
        T[-1,:] = T_top

	
        print("Iteration number :%s  Residuals U:%s T:%s" %(n,resU,resT))
	
        
    return u, v, p, T
