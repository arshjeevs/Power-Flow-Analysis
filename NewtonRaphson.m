clc;
clear;
format long;

% Load data (assumes data14.m script defines bus_data, line_data, shunt_data, tap_data)
data14;

Sbase = 100;
Vbase = 1;
bus_number = bus_data(:,1);
V_mag = bus_data(:,2);
V_angle = bus_data(:,3);
P_gen = bus_data(:,4) / Sbase;
Q_gen = bus_data(:,5) / Sbase;
P_load = bus_data(:,6) / Sbase;
Q_load = bus_data(:,7) / Sbase;
Q_min = bus_data(:,8) / Sbase;
Q_max = bus_data(:,9) / Sbase;

nbus = length(bus_number);
Ybus = Y_admi(line_data, shunt_data, tap_data, nbus);

% Define bus types
bus_type = ones(nbus,1) * 3; % All PQ initially
bus_type(1) = 1; % Slack
bus_type(2) = 2; % PV
bus_type(3) = 2; % PV

V = V_mag .* exp(1j * V_angle * pi/180);
P_spec = P_gen - P_load;
Q_spec = Q_gen - Q_load;

max_iter = 1000;
tol = 1e-6;

for iter = 1:max_iter
    S_calc = V .* conj(Ybus * V);
    P_calc = real(S_calc);
    Q_calc = imag(S_calc);

    % Enforce Q limits on PV buses
    for idx = 1:nbus
        if bus_type(idx) == 2
            if Q_calc(idx) < Q_min(idx)
                Q_spec(idx) = Q_min(idx);
                bus_type(idx) = 3;
            elseif Q_calc(idx) > Q_max(idx)
                Q_spec(idx) = Q_max(idx);
                bus_type(idx) = 3;
            end
        end
    end

    slack_bus = find(bus_type == 1);
    pv_buses = find(bus_type == 2);
    pq_buses = find(bus_type == 3);

    dP = P_spec - P_calc;
    dQ = Q_spec - Q_calc;

    mismatch = [dP([pv_buses; pq_buses]); dQ(pq_buses)];
    if max(abs(mismatch)) < tol
        fprintf('Converged in %d iterations.\n', iter);
        break;
    end

    J = calculate_jacobian(Ybus, V, bus_data, slack_bus, pv_buses, pq_buses);
    dx = J \ mismatch;

    dTheta = dx(1:length([pv_buses; pq_buses]));
    dV = dx(length(dTheta)+1:end);

    V_angle([pv_buses; pq_buses]) = V_angle([pv_buses; pq_buses]) + dTheta * 180/pi;
    V_mag(pq_buses) = V_mag(pq_buses) + dV;
    V = V_mag .* exp(1j * V_angle * pi/180);

    fprintf('Iteration: %d \n', iter);
end

if iter == max_iter
    fprintf('Did not converge in %d iterations.\n', max_iter);
end

fprintf('\nBus Voltage Results:\n');
fprintf('Bus #   Voltage (pu)   Angle (deg)  \n');
for i = 1:nbus
    fprintf('%4d %12.4f %12.4f \n', bus_number(i), abs(V(i)), angle(V(i))*180/pi);
end

fprintf('\nLine Power Flows:\n');
fprintf('From Bus   To Bus   P Flow (MW)   Q Flow (MVAR)\n');
for k = 1:size(line_data,1)
    from_bus = line_data(k,2);
    to_bus = line_data(k,3);
    I_from = (V(from_bus) - V(to_bus))*(-Ybus(from_bus,to_bus));
    S_from = V(from_bus) * conj(I_from);

    fprintf('%4d %9d %12.2f %12.2f\n', from_bus, to_bus, real(S_from)*Sbase, imag(S_from)*Sbase);
end

S = zeros(14, 14);
I_line = zeros(14, 14);
for i = 1:14
    for j = 1:14
        if i ~= j && Ybus(i, j) ~= 0
            I_line(i, j) = -Ybus(i, j) * (V(i) - V(j));
            S(i, j) = V(i) * conj(I_line(i, j));
        end
    end
end

S_gen = V(1) * conj(sum(I_line(1,:)));
disp('Power Generated at Slack Bus:');
disp(S_gen);

S_loss = sum(S(:));  
disp('Total Transmission Losses:');
disp(S_loss);