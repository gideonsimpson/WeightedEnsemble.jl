# routines for assembling a coarse model

"""
`build_coarse_transition_matrix`: Contruct a transition matrix amongst the bins (serial version).

### Arguments
* `mutation!` - an in place mutation function
* `bin_id` - bin identification function
* `x0_vals` - an array of starting values
* `bin0_vals` - an array of the bins corresponding to `x0_vals`
* `n_bins` - total number of bins
* `n_samples` - number of trials for each sample
"""
function build_coarse_transition_matrix(mutation!, bin_id, x0_vals, bin0_vals, n_bins, n_samples)

   X = similar(x0_vals[1]);
   K = spzeros(n_bins, n_bins);

   for k in 1:n_samples, l in 1:length(x0_vals)
      X .= deepcopy(x0_vals[l]);
      i = bin0_vals[l];
      mutation!(X);
      j = bin_id(X);
      K[i,j] +=1.0;
   end

   @. K = K/n_samples;

   return K
end


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
* `n_samples` - number of trials for each sample
"""

function pbuild_coarse_transition_matrix(mutation!, bin_id, x0_vals, bin0_vals, n_bins, n_samples)
   n_x0 = length(x0_vals);
   row_vals = SharedArray{Float64}(n_x0*n_samples);
   col_vals = SharedArray{Float64}(n_x0*n_samples);
   X = similar(x0_vals[1]);

   @sync @distributed for k in 1:n_samples
      for l in 1:n_x0
         X .= copy(x0_vals[l]);
         row_vals[n_x0*(k-1) + l] = bin0_vals[l];
         mutation!(X);
         col_vals[n_x0*(k-1) + l] = bin_id(X);
      end
   end

   K = sparse(row_vals,col_vals,ones(size(row_vals)),n_bins, n_bins);
   @. K = K/n_samples;
   return K
end

"""
`build_coarse_vectors`: Assemble the conditional expectation and 1- step
variance approximations one a coarser model, given the transition matrix, `K̃`,
and a coarse scale QoI function, `f̃`.

### Arguments
* `n_we_steps` - number of WE steps
* `K̃` - coarse scale transition matrix
* `f̃` - quantity of interest vector on the bin space
"""
function build_coarse_vectors(n_we_steps, K̃, f̃)
   n_bins = length(f̃);
   ṽ²_vals = [zeros(n_bins) for j in 1:n_we_steps];
   h̃_vals = [zeros(n_bins) for j in 1:n_we_steps+1];
   h̃ = deepcopy(f̃);
   K̃h̃ =similar(f̃);
   # vector of 1's
   e = ones(n_bins);
   @. h̃_vals[end] = h̃;

   # Use t+1 in indexing since Julia arrays start at 1
   for t in n_we_steps-1:-1:0

      K̃h̃ .= K̃ * h̃;
      # compute variance row by row
      for p in 1:n_bins
           ṽ²_vals[t+1][p] = K̃[p,:]⋅((h̃ .- e * (K̃h̃[p])).^2)
      end
      # update for next iterate
      @. h̃ = K̃h̃;
      @. h̃_vals[t+1] = h̃;
   end

   return h̃_vals, ṽ²_vals
end
