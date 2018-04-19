# Dynamic Time Warping
## The Problem
Assume having a reference signal ${y}_{R} = f \left( t \right)$.  
Given signals $\left \{ {y}_{i} = f \left( {s}_{i} \left( t \right) t \right ) \right \}$, namely temporally scaled versions of the same signal.  
How could on measure the similarity between the signals with no sensitivity to the "Tempo"?

***Example***

$$ \begin{align*}
{y}_{R} & = \sin \left( 2 \pi t \right), \; t \in [0, 1] \\
{y}_{1} & = \sin \left( 2 \pi 2 t \right), \; t \in [0, 0.5]
\end{align*} $$

![Figure 001](http://i.imgur.com/6gy1WOe.png)


**Intuition**
Think of a curve in 1D.  
Assume different people walking the path on different pace while the output is the height of ground the walk on (Or their head assuming they have the same height).



## Solution
> In Time Series Analysis, Dynamic Time Warping (DTW) is an algorithm for measuring similarity between two temporal sequences which may vary in speed.

Namely, the DTW provide us a way to measure similarity between signals with no sensitivity to their Time Scale.

### Where Can it Be Used?
 * Speech Recognition.
 * Audio Signature.
 * Shape Matching (Font Recognition).

### How Can It Be Done?
The answer is by creating a **Distance Matrix** and looking for the Optimal Match.

![Figure 002](http://i.imgur.com/G1hLSXr.png)

```matlab
vTimeSupport1 = [0:0.1:1].';
vTimeSupport2 = [0:0.2:1.6].';

vY1 = sin(2 * pi * vTimeSupport1);
vY2 = sin(2 * pi * vTimeSupport2);

mDistMtx = bsxfun(@minus, vY1, vY2.') .^ 2;
```


In the case above the Distance Matrix is given by $D \left( i, j \right) = { \left( {y}_{1} \left( i \right) - {y}_{2} \left( j \right) \right) }^{2}$.  
Yet one could use any Distance Measure.

### Optimal Match
Given the Distance Matrix the DTW is given by enforcing prior knowledge and extracting optimal path.
Paths are called **Warping Path** and defined as a sequence of coordinates $p = \left( {p}_{1}, {p}_{2}, \cdots, {p}_{L} \right)$.
The path cost is given by $pathCost = D \left( {p}_{1} \right) + D \left( {p}_{2} \right) + \cdots + D \left( {p}_{L} \right)$.
The Optimal (Minimum Cost Path) problem (Generally) can be solved by iterating all the possible paths. 

#### Classic Dynamic Time Warping
The **Warping Path** $ p $ satisfying:
 1. Boundary Condition - ${p}_{1} = \left( 1, 1 \right)$ and ${p}_{L} = \left( M, N \right)$.
 2. Monotonicity Condition - ${n}_{1} \leq {n}_{2} \leq \cdots \leq {n}_{L}$ and ${m}_{1} \leq {m}_{2} \leq \cdots \leq {m}_{L}$.
 3. Step Size Condition - ${p}_{l + 1} - {p}_{l} \in \left\{\left( 1, 0 \right), \left( 0, 1 \right), \left( 1, 1 \right) \right\}$.


Intuitively, assuming the signals start and end at the same time.  
Using the Walking Man metaphor, they start at the same point and end at the same point and each man must stand or walk forward, never go back.

The Optimal (Minimum Cost Path) problem under those constraints is a Dynamic Programming Problem given by:

$$ \gamma \left( i, j \right) = D \left( i, j \right ) + \begin{cases}
D \left( 1, 1 \right ) + D \left( 1, 2 \right ) + \cdots D \left( 1, j - 1 \right ) & \text{ if } i = 1 \\ 
D \left( 1, 1 \right ) + D \left( 2, 1 \right ) + \cdots D \left( i - 1, 1 \right ) & \text{ if } j = 1 \\ 
\min \left\{ \gamma \left( i - 1, j \right), \gamma \left( i, j - 1 \right), \gamma \left( i - 1, j - 1 \right) \right\} & \text{ if } i > 1, j > 1 
\end{cases} $$

Where $ \gamma \left( i, j \right) $ stands for the minimum cost to reach to point $ \left( i, j \right) $.

![Figure 003](http://i.imgur.com/YmaAl3h.png)

![Figure 004](http://i.imgur.com/fQs97t7.png)


#### Different Flavours
* Allowing shifted start / end point (Removing Constrain #001).  
  Update the "Initial Conditions" in the algorithm.
* Enforcing maximum distance between samples.
* Favoring Horizontal / Vertical / Diagonal Movement.
  For cases with prior knowledge on the path.
  Add weights for the "Base Points" for the current location.
* Different Step Size (Removing Constrain #003).
  Add "Base Points" to get to the current location.
* Multi Scale approach
  To find similar Sub Sequences.
* 

## Example

![Figure 005](http://i.imgur.com/tZOndT7.png)

```matlab
numSamples  = 60;
noiseStd    = 0.05;
shiftVal    = 5;

vX                  = [-1:0.05:1];
vGaussianDerivative = CalcGaussianGradient(3, 4);
vSigmoidSignal      = 1 ./ (1 + exp(-5 * vX));
vBaseSignal         = [zeros([1, 15]), vSigmoidSignal, ...
    ones([1, 20]), (1 + vGaussianDerivative), ones([1, 20])];
vBaseSignal         = 10 * vBaseSignal(:);
numSamplesBase      = size(vBaseSignal, 1);
vXBase              = [1:numSamplesBase].';

vY1Idx = sort(randperm(numSamplesBase, numSamples));
vY2Idx = sort(randperm(numSamplesBase, numSamples));

vY1 = vBaseSignal(vY1Idx);
vY2 = vBaseSignal(vY2Idx);

vY1 = vY1 + (noiseStd * randn([numSamples, 1]));
vY2 = vY2 + (noiseStd * randn([numSamples, 1]));
```

![Figure 006](http://i.imgur.com/IC9LkCV.png)

![Figure 007](http://i.imgur.com/3hAg0GM.png)

