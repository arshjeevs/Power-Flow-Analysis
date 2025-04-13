function J = calculate_jacobian(Ybus, V, bus_data, slack_bus, pv_buses, pq_buses)
    nbus = length(V);
    V_mag = abs(V);
    V_ang = angle(V);

    update_buses = sort([pv_buses; pq_buses]);
    n_update = length(update_buses);
    n_pq = length(pq_buses);

    J11 = zeros(n_update, n_update);
    J12 = zeros(n_update, n_pq);
    J21 = zeros(n_pq, n_update);
    J22 = zeros(n_pq, n_pq);

    S_calc = V .* conj(Ybus * V);
    P_calc = real(S_calc);
    Q_calc = imag(S_calc);

    for i_idx = 1:n_update
        i = update_buses(i_idx);
        for k_idx = 1:n_update
            k = update_buses(k_idx);
            if i == k
                J11(i_idx, k_idx) = -Q_calc(i) - imag(Ybus(i, i)) * V_mag(i)^2;
            else
                J11(i_idx, k_idx) = V_mag(i) * V_mag(k) * ...
                    (real(Ybus(i, k)) * sin(V_ang(i) - V_ang(k)) - ...
                     imag(Ybus(i, k)) * cos(V_ang(i) - V_ang(k)));
            end
        end
        for k_idx = 1:n_pq
            k = pq_buses(k_idx);
            if i == k
                J12(i_idx, k_idx) = P_calc(i) / V_mag(i) + real(Ybus(i, i)) * V_mag(i);
            else
                J12(i_idx, k_idx) = V_mag(i) * ...
                    (real(Ybus(i, k)) * cos(V_ang(i) - V_ang(k)) + ...
                     imag(Ybus(i, k)) * sin(V_ang(i) - V_ang(k)));
            end
        end
    end

    for i_idx = 1:n_pq
        i = pq_buses(i_idx);
        for k_idx = 1:n_update
            k = update_buses(k_idx);
            if i == k
                J21(i_idx, k_idx) = P_calc(i) - real(Ybus(i, i)) * V_mag(i)^2;
            else
                J21(i_idx, k_idx) = -V_mag(i) * V_mag(k) * ...
                    (real(Ybus(i, k)) * cos(V_ang(i) - V_ang(k)) + ...
                     imag(Ybus(i, k)) * sin(V_ang(i) - V_ang(k)));
            end
        end
        for k_idx = 1:n_pq
            k = pq_buses(k_idx);
            if i == k
                J22(i_idx, k_idx) = Q_calc(i) / V_mag(i) - imag(Ybus(i, i)) * V_mag(i);
            else
                J22(i_idx, k_idx) = V_mag(i) * ...
                    (real(Ybus(i, k)) * sin(V_ang(i) - V_ang(k)) - ...
                     imag(Ybus(i, k)) * cos(V_ang(i) - V_ang(k)));
            end
        end
    end

    J = [J11 J12; J21 J22];
end
