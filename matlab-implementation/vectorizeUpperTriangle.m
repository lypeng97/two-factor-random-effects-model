function output = vectorizeUpperTriangle(inputMatrix)
    % Check the dimensions of inputMatrix
    dims = size(inputMatrix);
    
    if length(dims) == 3
        % If input is 3D: subject x region x region
        [num_subjects, num_regions, ~] = size(inputMatrix);
        
        % Calculate the number of upper triangle elements
        num_upper_triangle_values = num_regions * (num_regions - 1) / 2;

        % Preallocate a matrix to hold the upper triangle values for each subject
        output = zeros(num_subjects, num_upper_triangle_values);

        % Loop through each subject and extract the upper triangle
        for subject = 1:num_subjects
            % Extract the subject's 2D connectivity matrix
            subject_connectivity_strength = squeeze(inputMatrix(subject, :, :));

            % Extract the upper triangle (excluding diagonal elements)
            upper_triangle_indices = triu(true(size(subject_connectivity_strength)), 1);
            upper_triangle_values = subject_connectivity_strength(upper_triangle_indices);

            % Store the upper triangle values in the corresponding row of the matrix
            output(subject, :) = upper_triangle_values;
        end

    elseif length(dims) == 2 && dims(1) == dims(2)
        % If input is a single 2D symmetric matrix
        num_regions = dims(1);
        
        % Calculate the number of upper triangle elements
        num_upper_triangle_values = num_regions * (num_regions - 1) / 2;
        
        % Extract the upper triangle (excluding diagonal elements)
        upper_triangle_indices = triu(true(num_regions), 1);
        output = inputMatrix(upper_triangle_indices)';

    else
        error('Input must be a 3D matrix (subject x region x region) or a 2D symmetric matrix.');
    end
end
