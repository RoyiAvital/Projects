[![Visitors](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FRoyiAvital%2FStackExchangeCodes&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=Visitors+%28Daily+%2F+Total%29&edge_flat=false)](https://github.com/RoyiAvital/Julia100Exercises)
[![Royi Avital](https://i.imgur.com/ghq7NUE.png)](https://github.com/RoyiAvital/StackExchangeCodes)

# Linear Segmentation
This document describe a [_Dynamic Programming_](https://en.wikipedia.org/wiki/Dynamic_programming) based method to solve the _linear segmentation_ problem.

> Notebook by:  
> - Royi Avital RoyiAvital@yahoo.com

> Attribution
> 

References:

 1. [RcppDynProg](https://github.com/WinVector/RcppDynProg)  
    Defined the problem. They don't define their solution (There is code I haven't looked at).

Remarks:

 1. Inspired by [`LinearSegmentation.jl`](https://github.com/stelmo/LinearSegmentation.jl).

To Do:

 1. Add support for affinity (Score) in addition to distance.

## Revision History

| Version | Date       | User        |Content / Changes                                                                         |
|---------|------------|-------------|------------------------------------------------------------------------------------------|
| 0.1.000 | 12/07/2023 | Royi Avital | First version                                                                            |

## Problem Description

Given a regression data set $\left\{ {x}_{i}, {y}_{i} \right\}_{i = 1}^{N}$ find the optimal piece wise linear model for the data.  

![Piece Wise Linear Model](https://imgur.com/QKa8mZX.png)


One could easily estimate the linear function parameters given the different segments.  
The challenge is to infer the segments without any side / prior information.  

## Proposed Solution

The idea of the solution is as following:

 * Define the distance between data points $d \left( \left( {x}_{i}, {y}_{i} \right), \left( {x}_{j}, {y}_{j} \right) \right)$.
 * Define the cost of a segment $s \left[ i, j \right]$ to be the distance between its edge points $\operatorname{cost} \left( s \left[ i, j \right] \right) = d \left( \left( {x}_{i}, {y}_{i} \right), \left( {x}_{j}, {y}_{j} \right) \right)$.
 * Partition the data into non overlapping segments which minimizes the total cost of the segments.

The problem above is a specific case of the more general partitioning problem (In general a _NP Hard_ problem).  

### Dynamic Programming Solution

In order to solve the problem, the data is modeled as a cost matrix $\boldsymbol{C}$:

$$ \boldsymbol{C}_{i, j} = \begin{cases} d \left( i, j \right ) & \text{ if } i \leq j \\ \infty & \text{ if } i > j \end{cases} $$

Where $d \left( i, j \right ) = d \left( \left( {x}_{i}, {y}_{i} \right), \left( {x}_{j}, {y}_{j} \right) \right)$.

For example, assume the following data: $\boldsymbol{d} = {\left[ 1, 1.1, 0.9, 7, 8, 7.5, 4, 3.6, 4.4 \right]}^{T}$ with $d \left( i, j \right) = {\left( {d}_{i} - {d}_{j} \right)}^{2}$.  
Then the cost matrix $\boldsymbol{C}$ is given by:

![Cost Function](https://i.imgur.com/P2jpZEW.png)

One way to solve it is by defining the matrix which defines the _Dynamic Programming_ path:

$$ \boldsymbol{S}_{i, j} = \begin{cases} \min_{k \in \left\{ 1, 2, \ldots, i - 1 \right\}} \boldsymbol{S}_{k, i - 1} + \boldsymbol{C}_{i, j} & \text{ if } i < j \\ \infty & \text{ if } i => j \end{cases} $$

Where:

 * $\boldsymbol{S}_{i, j}$ is the cost of having a segment $i \to j$ with the previous segment finishes at $i - 1$.
 * From all previous segments that finish at $i - 1$ the segment is appended to the one the minimum cost.
 * The recursive formulation suggests all combinations are evaluated to the optimal solution.
 * Since the $i$ -th segment must starts at least from the point $i$ means the solution is optimal.
 * In addition, to recover the actual segments one need to keep the $ k $ parameter per segment.
 * The optimal path of segments starts at lowest cost path which at the last column.  
   Namely at $\boldsymbol{S}_{i, N}, \; i = \arg \min_{k} \boldsymbol{S}_{i, N}$ where $N$ is the number of points in data.


|![Segment Matrix](https://i.imgur.com/p3G9nlP.png)|![Path Matrix](https://i.imgur.com/KERnSSd.png)|
|--------------------------------------------------|-----------------------------------------------|

The functions to calculate the 

```julia
function SolveMinCostPartitionIntervals( mD :: Matrix{T}, maxPartitions :: S  ) where {T, S <: Integer}

    numSamples = size(mD, 1);
    maxPartitions = min(maxPartitions, numSamples);
    mS = maximum(mD) .* ones(T, maxPartitions, numSamples); #<! Cost per segment
    mP = zeros(Int, maxPartitions, numSamples); #<! Path matrix

    mS[1, :] .= mD[1, :];
    mP[1, :] .= 1;
    for ii ∈ 2:maxPartitions
        for jj ∈ (ii + 1):numSamples
            minCost = 1e6;
            kkMin = 1;
            for kk ∈ 1:(ii - 1)
                currCost = mS[kk, ii - 1] + mD[ii, jj];
                if (currCost < minCost)
                    kkMin   = kk;
                    minCost = currCost;
                end
            mS[ii, jj] = minCost;
            mP[ii, jj] = kkMin;
            end
        end
    end

    return mS, mP;

end

function ExtractPath( mS :: Matrix{T}, mP :: Matrix{S} ) where {T, S <: Integer}

    numRows, numCols    = size(mS);
    vS                  = Vector{Vector{Int}}(undef, 0);

    startIdx = argmin(mS[:, numCols]);
    endIdx   = numCols;

    while ((startIdx > 0) && (endIdx > 0))
        prepend!(vS, [[startIdx, endIdx]]);
        colIdx      = mP[startIdx, endIdx];
        endIdx      = startIdx - 1;
        startIdx    = colIdx;
    end

    return vS;

end
```

With the output being colored:

![Segmented Data](https://i.imgur.com/DuQjbv7.png)

### Linear Segmentation

To utilize the method in the context of _piece wise linear regression_ the following steps are taken:

 * Define the _Distance Matrix_ (_Cost Matrix_) $\boldsymbol{D}$ as the MSE score for the segment $\left[ i, j \right]$.  
   Namely, all samples within the segment are taken to build a local linear function and the MSE between the function and the samples is used as the cost.  
   In the implementation of such distance function one can set the length of the minimal and maximal segments.  
   Specifically to take care of the case of 1 and 2 points.
 * The distance matrix is analyzed for the best segment partitioning as above.
 * The output is a vector of 2D vectors which define the segments.


> **Note**  
> The distance function is implemented to support any loss function.  
> For such use case, the [R2 Score](https://en.wikipedia.org/wiki/Coefficient_of_determination), being normalized, make sense.  
> To generate a loss function out of it requires a simple transformation.

### Example

In this example a signal is generated as an harmonic signal with time variant phase and frequency added to a piece wise linear data.

![Piece Wise Linear Data with Harmonic Noise](https://i.imgur.com/45OAMwV.png)

With the selection of `minSegLen = 10.0` the segmentation resulted in:

![Segmented Data](https://i.imgur.com/prQJ2n0.png)

## Replicate the Figures

 * Activate the Julia environment.
 * Run `LinearSegmentationData.jl` to generate data.
 * Run `LinearSegmentation.jl` to generate figures.
 * The functions are in `LinearSegmentationFun.jl`.



