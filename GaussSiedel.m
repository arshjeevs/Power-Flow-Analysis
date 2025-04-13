clc; clear; close all;
format long;
data14;
nbus=length(bus_data(:,1));
function [V, Q, iter] = gauss_seidel_power_flow(Ybus, V, P, Q,Q_max,Q_min, PQ_buses, PV_buses, tol, max_iter)
    n = length(V); 
    iter = 0;

    while iter < max_iter
        V_old = V; 
        iter = iter + 1;
        for k = 1:length(PQ_buses)
            i = PQ_buses(k);
            sumYV = 0;
            for j = 1:n
                if j ~= i
                    sumYV = sumYV + Ybus(i, j) * V(j);
                end
            end
            S_i = P(i) - 1j * Q(i);
            V(i) = (1 / Ybus(i, i)) * ((S_i / conj(V_old(i))) - sumYV);
            
        end
        %fprintf('Iteration: %d \n', iter);
        %disp(V)
        for k = 1:length(PV_buses)
            i = PV_buses(k);
            sumYV = 0;
            for j = 1:n
                if j ~= i
                    sumYV = sumYV + Ybus(i, j) * V(j);
                end
            end
            
            Q_calc(i) = -imag(conj(V_old(i)) * (Ybus(i, i) * V_old(i) + sumYV));
            if Q_calc(i) < Q_min(i)
                    Q_calc(i) = Q_min(i);
                    PQ_buses(end+1)=i;
                    PV_buses(k)=[];% Switch to PQ bus if limit violated

            elseif Q_calc(i) > Q_max(i)
                    Q_calc(i) = Q_max(i);
                    PQ_buses(end+1)=i;
                    PV_buses(k)=[]; % Switch to PQ bus if limit violated
            end
            Q(i)=Q_calc(i);
            S_i = P(i) - 1j * Q(i);
            V_new = (1 / Ybus(i, i)) * ((S_i / conj(V_old(i))) - sumYV);
            V(i) = abs(V_old(i)) * (V_new / abs(V_new));           
        end

        if max(abs(V - V_old)) < tol
            fprintf('Power Flow Converged in %d iterations.\n', iter);
            return;
        end
    end
    fprintf('Did not converge within %d iterations.\n', max_iter);
end


%%%%%%%%%%    
Sbase = 100;
bus_number = bus_data(:,1);
V_mag = bus_data(:,2);
V_angle = bus_data(:,3);
P_gen = bus_data(:,4) / Sbase;
Q_gen = bus_data(:,5) / Sbase;
P_load = bus_data(:,6) / Sbase;
Q_load = bus_data(:,7) / Sbase;
Q_min = bus_data(:,8) / Sbase;
Q_max = bus_data(:,9) / Sbase;


Ybus = Y_admi(line_data,shunt_data, tap_data, nbus);  

P= P_gen-P_load;     
Q= Q_gen- Q_load;  

PQ_buses = [4,5,6,7,8,9,10,11,12,13,14]; 
PV_buses = [2,3];    

tol = 1e-6;
max_iter = 1000;

[V_final, Q_final, iterations] = gauss_seidel_power_flow(Ybus, V_mag, P, Q,Q_max,Q_min, PQ_buses, PV_buses, tol, max_iter);

disp('Final Bus Voltages:');
disp(V_final);

S = zeros(14, 14);
I_line = zeros(14, 14);
for i = 1:14
    for j = 1:14
        if i ~= j && Ybus(i, j) ~= 0
            I_line(i, j) = -Ybus(i, j) * (V_final(i) - V_final(j));  % Line current
            S(i, j) = V_final(i) * conj(I_line(i, j));  % Power Flow
        end
    end
end

fprintf('\nBus Voltage Results:\n');
fprintf('Bus #   Voltage (pu)   Angle (deg)  \n');
for i = 1:nbus
    fprintf('%4d %12.4f %12.4f \n', bus_number(i), abs(V_final(i)), angle(V_final(i))*180/pi);
end

    fprintf('\nLine Power Flows:\n');
    fprintf('From Bus   To Bus   P Flow (MW)   Q Flow (MVAR)\n');
    for k = 1:size(line_data,1)
        from_bus = line_data(k,2);
        to_bus = line_data(k,3);
        R = line_data(k,4);
        X = line_data(k,5);
        B = line_data(k,6);
        
        % Check if this line is a transformer
        is_transformer = 0;
        tap = 1;
        for m = 1:size(tap_data,1)
            if (from_bus == tap_data(m,1) && to_bus == tap_data(m,2)) || ...
               (from_bus == tap_data(m,2) && to_bus == tap_data(m,1))
                tap = tap_data(m,3);
                is_transformer = 1;
                break;
            end
        end
        
        if is_transformer
            if from_bus == tap_data(m,1) % From is primary side
                I_from = (V_final(from_bus)/tap - V_final(to_bus)) / (R + 1j*X);
                S_from = V_final(from_bus)/tap * conj(I_from);
            else % From is secondary side
                I_from = (V_final(from_bus) - V_final(to_bus)/tap) / (R + 1j*X);
                S_from = V_final(from_bus) * conj(I_from);
            end
            I_to = -I_from;
            S_to = V_final(to_bus) * conj(I_to);
        else
            I_from = (V_final(from_bus) - V_final(to_bus)) / (R + 1j*X) + V_final(from_bus) * (1j*B/2);
            S_from = V_final(from_bus) * conj(I_from);
            I_to = (V_final(to_bus) - V_final(from_bus)) / (R + 1j*X) + V_final(to_bus) * (1j*B/2);
            S_to = V_final(to_bus) * conj(I_to);
        end
        
        fprintf('%4d %9d %12.2f %12.2f\n', from_bus, to_bus, real(S_from)*Sbase, imag(S_from)*Sbase);
    end





S_gen = V_final(1) * conj(sum(I_line(1,:)));
disp('Power flowing out of Slack Bus:');
disp(S_gen);

S_loss = sum(S(:));  
disp('Total Transmission Losses:');
disp(S_loss);
