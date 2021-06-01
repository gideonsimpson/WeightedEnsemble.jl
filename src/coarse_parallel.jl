"""
`pbuild_coarse_transition_matrix`: Contruct a transition matrix amongst the
bins in parallel.  It assumed that a pool of workers has already been
contructed.  Returns a sparse matrix.

### Arguments
* `mutation!` - an in place mutation function
* `bin_id` - bin identification function
* `x0_vals` - an array of starting values
* `bin0_vals` - an array of the bins corresponding to `x0_vals`
* `n_bins` - total number of bins
* `n_trials` - number of independent trials for each x0 starting value
"""

function pbuild_coarse_transition_matrix(mutation!, bin_id, x0_vals, n_bins, n_trials)
   n_x0 = length(x0_vals);
   row_vals = SharedArray{Int}(n_x0*n_trials);
   col_vals = SharedArray{Int}(n_x0*n_trials);
   X = similar(x0_vals[1]);

   @sync @distributed for k in 1:n_trials
      for l in 1:n_x0
         X .= deepcopy(x0_vals[l]);
         row_vals[n_x0*(k-1) + l] = bin_id(X);
         mutation!(X);
         col_vals[n_x0*(k-1) + l] = bin_id(X);
      end
   end

   K = sparse(row_vals,col_vals,ones(size(row_vals)),n_bins, n_bins);
   for i in 1:n_bins
      rowZ = sum(K[i,:]);
      @. K[i,:] = K[i,:]/rowZ; 
   end
   return K
end
