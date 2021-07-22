# Goal
This is an example implementation of the transmission process using polar 
coding. It is intended to complement the original paper on polar coding [1].

# Remarks
main_transmit_on_BEC.m is the main script.
It simulates the transmission process on a BEC(EPSILON) at the given RATE 
and with the given BLOCKLENGTH. 

There are two versions of the successive cancellation decoder. Their 
output are exactly similar, however the 'naive' version is less efficient
than the optimized one (complexity O(N^2) vs O(NlogN)).

For questions or comments, please send a mail to fsabatier0@gmail.com

# References
[1] E. Arikan, Channel polarization: a method for constructing 
capacity-achieving codes for symmetric binary-input memoryless channels, 
IEEE Trans. Inf. Theory, vol. 55, no. 7, pp. 3051-3073, July 2009.


