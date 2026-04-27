function Y = extract_measurements(R)

fn = fieldnames(R);
Y = R.(fn{1});

if size(Y,1) ~= 6
    Y = Y';
end

end