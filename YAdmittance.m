function [Ybus] = Y_admi(line_data,shunt_data, tap_data, nbus)
%% Y_BUS Building
Ybus = zeros(nbus, nbus);
    
    % Build Ybus matrix
    for k = 1:size(line_data,1)
        from_bus = line_data(k,2);
        to_bus = line_data(k,3);
        R = line_data(k,4);
        X = line_data(k,5);
        B = line_data(k,6);
        
        % Check if this line is a transformer
        is_transformer = 0;
        for m = 1:size(tap_data,1)
            if (from_bus == tap_data(m,1) && to_bus == tap_data(m,2)) || ...
               (from_bus == tap_data(m,2) && to_bus == tap_data(m,1))
                tap = tap_data(m,3);
                is_transformer = 1;
                break;
            end
        end
        
        if is_transformer
            % Transformer admittance
            y = 1 / (R + 1j*X);
            Ybus(from_bus, to_bus) = Ybus(from_bus, to_bus) - y/tap;
            Ybus(to_bus, from_bus) = Ybus(to_bus, from_bus) - y/tap;
            Ybus(from_bus, from_bus) = Ybus(from_bus, from_bus) + y/(tap^2);
            Ybus(to_bus, to_bus) = Ybus(to_bus, to_bus) + y;
        else
            % Regular line admittance
            y = 1 / (R + 1j*X);
            Ybus(from_bus, to_bus) = Ybus(from_bus, to_bus) - y;
            Ybus(to_bus, from_bus) = Ybus(to_bus, from_bus) - y;
            Ybus(from_bus, from_bus) = Ybus(from_bus, from_bus) + y + 1j*B/2;
            Ybus(to_bus, to_bus) = Ybus(to_bus, to_bus) + y + 1j*B/2;
        end
    end
    
    % Add shunt capacitors
    for k = 1:size(shunt_data,1)
        bus = shunt_data(k,1);
        B_shunt = shunt_data(k,2);
        Ybus(bus, bus) = Ybus(bus, bus) + 1j*B_shunt;
    end

