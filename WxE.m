
tspan = [0,10];

# love equation function
function [t,x] = system_simulation(A, tspan, y0)
    ode_sys = @(t,x) [A(1,1)*x(1)+A(1,2)*x(2);A(2,1)*x(1)+A(2,2)*x(2)]; # Definition system
    [t,x] = ode23(ode_sys, tspan, y0); # Resolution system
endfunction

# Calculer les isoclines
function [isocline_1,isocline_2] = compute_isoclines(A,line_range)
    isocline_1 = -(A(1,1)/A(1,2)) * line_range;
    isocline_2 = -(A(2,1)/A(2,2)) * line_range;
endfunction

function [line_range,isocline_1,isocline_2] = plot_isoclines(A)
    line_range = -1.5:.1:1.5;
    [isocline_1,isocline_2] = compute_isoclines(A,line_range);
endfunction

# eigenlines
function [eigenline_1,eigenline_2,V] = compute_eigenlines(A,line_range)
    [V,L] = eig(A);
    eigenline_1 = (V(2,1)/V(1,1)) * line_range;
    eigenline_2 = (V(2,2)/V(1,2)) * line_range;
endfunction

function [line_range,eigenline_1,eigenline_2] = plot_eigenlines(A)
    line_range = -1.5:.1:1.5;
    [eigenline_1,eigenline_2,V] = compute_eigenlines(A,line_range);
endfunction

# 4. Portrait de phase complet
function [x1,x2,x1p,x2p] = plot_portrait_phase_complete(A,a,b)
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
    plot(x1range,isocline_1,"linewidth",1,'LineStyle','--');
    plot(x1range,isocline_2,"linewidth",1,'LineStyle','--');
    # Vecteurs propres = droites invariantes
    plot(x1range,eigenline_1,"linewidth",1);
    plot(x1range,eigenline_2,"linewidth",1);
    titl = ["Phase portrait for a=",num2str(a)," b=",num2str(b)];
    title(titl)
    legend("field","isocline_1","isocline_2","v_1","v_2","location","south","orientation", "horizontal");
endfunction

function system_simulation_plot(A,tspan,y0,figure_number,a,b)
    [t,x] = system_simulation(A,tspan,y0);
    figure(figure_number);
    plot(t,x)
    titl = ["Evolution of w(t) and e(t) for t \in [0, 10], a=",num2str(a),", b=",num2str(b)];
    title(titl)
    xlabel("t")
    ylabel("x")
    legend("w(t)","e(t)","location","south","orientation", "horizontal");
endfunction

%function system_simulation_and_portrait_phase_and_plot(tspan,coefficient_matrices,y0)
%function system_simulation_and_portrait_phase_and_plot(tspan,y0)
function system_simulation_and_portrait_phase_and_plot(tspan)
    disp('system_simulation_and_portrait_phase_and_plot()')
    figure_number = 11
    initial_conditions_list = {[1,1]
                              %,[-1,1]
                              %,[-2,2]
                              }
    listOfTuples = {[-0.15,0.9]
                    %,[-0.15,-0.9]
                    %,[0.15,0.9]
                    %,[0.15,-0.9]
                    %,[1,0.1]
                    %,[0.1,0.1]
                    %,[-0.1,-0.1]
                    }

    fprintf('listOfTuples length : %d .\n',length(listOfTuples));
    fprintf('initial_conditions_list length : %d .\n',length(initial_conditions_list));

    for coefficient_matrix_index = 1:4
        for i = 1:length(listOfTuples)
            for ic = 1:length(initial_conditions_list)
                % Get the current initial_conditions
                fprintf('initial_conditions : %d .\n',initial_conditions_list{ic});

                y0 = initial_conditions_list{ic};
                %y0 = initial_conditions_list(ic);
                % Get the current a & b
                currentTuple = listOfTuples{i};
                %currentTuple = listOfTuples(i);

                % Extract the values from the tuple into two variables
                a = currentTuple(1);
                b = currentTuple(2);

                % print the values
                %fprintf('a : %d .\n', a);
                %fprintf('b : %d .\n', b);

                if coefficient_matrix_index == 1
                    #1
                    A = [a b;
                          0 0];
                elseif coefficient_matrix_index == 2
                    #2
                    A = [a b;
                          b a];
                elseif coefficient_matrix_index == 3
                    #3
                    A = [a b;
                          -b -a];
                elseif coefficient_matrix_index == 4
                    #4
                    A = [a 0;
                          b 0];
                end
                disp(A)

                system_simulation_plot(A,tspan,y0,figure_number,a,b);
                figure_number = figure_number + 1;
                # display portrait de phases
                %figure figure_number;
                figure(figure_number);
                [x1,x2,x1p,x2p] = plot_portrait_phase_complete(A,a,b);
                figure_number = figure_number + 1;
            endfor
        endfor
    endfor
endfunction


%system_simulation_and_portrait_phase_and_plot(tspan,coefficient_matrices,y0)
%system_simulation_and_portrait_phase_and_plot(tspan,y0)
system_simulation_and_portrait_phase_and_plot(tspan)
