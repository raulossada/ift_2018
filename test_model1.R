#####
library('deSolve');
library('reshape2');
library('ggplot2');

# time sequence 
time <- seq(from=0, to=300, by=0.01);

# parameters
parameters <- c(r_Z=1, K_z=1e10, alpha=0, g=1, m1=0.75, mu_Pz=1, s=1, i=1, 
                r_P=1, K_P=100, theta=1, m2=0.30, f2=0.70,
                r_Pz=1, K_Pz=100, f1=0.25, 
                delta=1, mu_D=1, 
                k1=1, k2=1, mu_C=1);

# initial conditions
state <- c(Z=1, P=1, Pz=0, D=1, C=1);

# model
model1 <- function(t, state, parameters){
  with( as.list( c(state, parameters) ), {
    
    dZ  <- r_Z*Z*(1-Z/K_z) + (1-alpha)*g*(m1*D + mu_Pz)*Pz - s*Z*P - i*Z*C;
    dP  <- r_P*P*(1-P/K_P) - theta*s*Z*P - (m2+f2)*P*D;
    dPz <- r_Pz*Pz*(1-Pz/K_Pz) + theta*s*Z*P - (m1+f1)*Pz*D;
    dD  <- delta*( (m2+f2)*P + (m1+f1)*Pz )*D - mu_D*D;
    dC  <- k1*i*Z*C + k2*alpha*g*(m1*D + mu_Pz)*Pz*C - mu_C*C;
    
    return( list( c(dZ, dP, dPz, dD, dC) ) );
  })
}

# Numerical solution
res <- ode(y=state, times=time, func=model1, parms=parameters, method='lsoda');



#####
# Plot
res_df <- as.data.frame(x=res);
res_m <- melt(data=res_df, id.vars='time');

###
p <- ggplot(data=res_m, mapping=aes(x=time, y=value, color=variable, colour=c('green', 'orange', 'blue', 'purple', 'red') ) ) + 
  geom_point(size=0.5 ) + 
  labs(y='Populations') + 
  theme_gdocs() + 
  scale_color_gdocs();
print(p);


###
Z_m <- res_m[which(res_m$variable=='Z'), ];
p_Z <- ggplot(data=Z_m, mapping=aes(x=time, y=value) ) + 
  geom_point(size=0.5, col='green') + 
  labs(y='Number of Zoochlorella (Z)');
print(p_Z);

###
P_m <- res_m[which(res_m$variable=='P'), ];
p_P <- ggplot(data=P_m, mapping=aes(x=time, y=value) ) + 
  geom_point(size=0.5, col='orange') + 
  labs(y='Number of Paramecium spp (P)');
print(p_P);

###
Pz_m <- res_m[which(res_m$variable=='Pz'), ];
p_Pz <- ggplot(data=Pz_m, mapping=aes(x=time, y=value) ) + 
  geom_point(size=0.5, col='blue') + 
  labs(y='Number of Paramecium spp with Zoochlorella (Pz)');
print(p_Pz);

###
D_m <- res_m[which(res_m$variable=='D'), ];
p_D <- ggplot(data=D_m, mapping=aes(x=time, y=value) ) + 
  geom_point(size=0.5, col='purple') + 
  labs(y='Number of Didinium nasutum (D)');
print(p_D);

###
C_m <- res_m[which(res_m$variable=='C'), ];
p_C <- ggplot(data=C_m, mapping=aes(x=time, y=value, color=variable) ) + 
  geom_point(size=0.5, col='red') + 
  labs(y='Number of Chlorovirus (C)');
print(p_C);
