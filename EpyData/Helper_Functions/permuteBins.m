function shuffled_vector = permuteBins(vector, bin_size)
    % This function divides a given vector into bins of size `bin_size`,
    % randomly shuffles these bins, and reconstructs the vector.
    % If the vector length is not a perfect multiple of `bin_size`,
    % the remaining elements are appended from a randomly chosen bin.

    % Determine the number of full bins that can be formed
    n_bins = floor(length(vector) / bin_size);

    % Reshape the vector into a matrix where each column is a bin
    reshaped_vector = reshape(vector(1:n_bins*bin_size), bin_size, n_bins);

    % Randomly shuffle the order of the bins
    permuted_bins = reshaped_vector(:, randperm(n_bins));

    % Flatten the shuffled matrix back into a vector
    shuffled_vector = permuted_bins(:);

    % Handle leftover elements
    n_leftover = length(vector) - length(shuffled_vector);
    
    if n_leftover > 0
        random_integer = randi(n_bins);
        leftover_data = reshaped_vector(1:n_leftover, random_integer);
        shuffled_vector = [shuffled_vector; leftover_data];
    end
end