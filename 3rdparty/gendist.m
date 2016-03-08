function T = gendist(P,N,M,varargin)
%
%GENDIST - generate random numbers according to a discrete probability
%distribution.
%Tristan Ursell, (c) 2012
%
%T = gendist(P,N,M)
%T = gendist(P,N,M,'plot')
%
%The function gendist(P,N,M) takes in a positive vector P whose values
%form a discrete probability distribution for the indices of P. The
%function outputs an N x M matrix of integers corresponding to the indices
%of P chosen at random from the given underlying distribution.
%
%P will be normalized, if it is not normalized already. Both N and M must
%be greater than or equal to 1.  The optional parameter 'plot' creates a
%plot displaying the input distribution in red and the generated points as
%a blue histogram.
%
%Conceptual EXAMPLE:
%
%If P = [0.2 0.4 0.4] (note sum(P)=1), then T can only take on values of 1,
%2 or 3, corresponding to the possible indices of P.  If one calls 
%gendist(P,1,10), then on average the output T should contain two 1's, four
%2's and four 3's, in accordance with the values of P, and a possible 
%output might look like:
%
%T = gendist(P,1,10)
%T = 2 2 2 3 3 3 1 3 1 3
%
%If, for example, P is a probability distribution for position, with
%positions X = [5 10 15] (does not have to be linearly spaced), then the
%set of generated positions is X(T).
%
%EXAMPLE 1:
%
%P = rand(1,50);
%T = gendist(P,100,1000,'plot');
%
%EXAMPLE 2:
%
%X = -3:0.1:3;
%P = 1+sin(X);
%T = gendist(P,100,1000,'plot');
%Xrand = X(T);
%
%Note:
%size(T) = 100 1000
%size(Xrand) = 100 1000
%

%Thanks to Derek O'Connor for tips on speeding up the code.

%check for input errors
if and(nargin~=3,nargin~=4)
    error('Error:  Invalid number of input arguments.')
end

if min(P)<0
    error('Error:  All elements of first argument, P, must be positive.')
end

if or(N<1,M<1)
    error('Error:  Output matrix dimensions must be greater than or equal to one.')
end

%normalize P
Pnorm=[0 P]/sum(P);

%create cumlative distribution
Pcum=cumsum(Pnorm);

%create random matrix
N=round(N);
M=round(M);
R=rand(1,N*M);

%calculate T output matrix
V=1:length(P);
[~,inds] = histc(R,Pcum); 
T = V(inds);

%shape into output matrix
T=reshape(T,N,M);

%if desired, output plot
if nargin==4
    if strcmp(varargin{1},'plot')
        Pfreq=N*M*P/sum(P);
        figure;
        hold on
        hist(T(T>0),1:length(P))
        plot(Pfreq,'r-o')
        ylabel('Frequency')
        xlabel('P-vector Index')
        axis tight
        box on
    end
end
