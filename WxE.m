
y0 = [-1,1];
tspan = [0,10];
a = -0.15;
b = 0.9;

# coefficient matrix
#1
A1 = [a b;
      0 0];
#2
A2 = [a b;
      b a];
#3
A3 = [a b;
      -b -a];
#4
A4 = [a 0;
      b 0];

# love equation function
function [t,x] = system_simulation(A, tspan, y0)
    ode_sys = @(t,x) [A(1,1)*x(1)+A(1,2)*x(2);A(2,1)*x(1)+A(2,2)*x(2)]; # Définition système
    [t,x] = ode23(ode_sys, tspan, y0); # Résolution système
endfunction

A = A4; # change Ax to plot another equation
[t,x] = system_simulation(A, tspan, y0);
figure 1
plot(t,x)


# Calculer les isoclines
function [isocline_1,isocline_2] = compute_isoclines(A,line_range)
    isocline_1 = -(A(1,1)/A(1,2)) * line_range;
    isocline_2 = -(A(2,1)/A(2,2)) * line_range;
endfunction

function [line_range,isocline_1,isocline_2] = plot_isoclines(A)
    line_range = -1.5:.1:1.5;
    [isocline_1,isocline_2] = compute_isoclines(A,line_range);
endfunction

[line_range,isocline_1,isocline_2] = plot_isoclines(A);
figure 2
plot(line_range,isocline_1,"linewidth",1);
hold on;
plot(line_range,isocline_2,"linewidth",1);
legend("isocline_1","isocline_2","location","south");

# 4. Portrait de phase
function [x1,x2,x1p,x2p] = plot_portrait_phase(A)
    #grid for plotting
    x1range=-1.5:.1:1.5;
    x2range=-1.5:.1:1.5;
    [x1,x2] = meshgrid(x1range, x2range);

    # Define the system to plot (based on matrix A)
    x1p = A(1,1)*x1+A(1,2)*x2;
    x2p = A(2,1)*x1+A(2,2)*x2;

    #Normalize values for plotting
    norms=sqrt(x1p.^2+x2p.^2);
    # Vector field plot
    hold on;
    quiver(x1,x2,x1p./norms,x2p./norms,0.5);
endfunction

[x1,x2,x1p,x2p] = plot_portrait_phase(A);

# 4. Portrait de phase complet
function [x1,x2,x1p,x2p] = plot_portrait_phase_complete(A)
    #grid for plotting
    x1range=-1.5:.1:1.5;
    x2range=-1.5:.1:1.5;
    [x1,x2] = meshgrid(x1range, x2range);

    # Define the system to plot (based on matrix A)
    x1p = A(1,1)*x1+A(1,2)*x2;
    x2p = A(2,1)*x1+A(2,2)*x2;

    #Normalize values for plotting
    norms=sqrt(x1p.^2+x2p.^2);

    [eigenline_1,eigenline_2,V] = compute_eigenlines(A,x1range);
    [isocline_1,isocline_2] = compute_isoclines(A,x1range);

    # Vector field plot
    hold on;
    quiver(x1,x2,x1p./norms,x2p./norms,0.5);
    # Isoclines
    plot(x1range,isocline_1,"linewidth",1);
    plot(x1range,isocline_2,"linewidth",1);
    # Vecteurs propres
    plot(x1range,eigenline_1,"linewidth",1);
    plot(x1range,eigenline_2,"linewidth",1);
    legend("field","v_1","v_2","isocline_1","isocline_2","location","south","orientation", "horizontal");
endfunction

[x1,x2,x1p,x2p] = plot_portrait_phase_complete(A);













