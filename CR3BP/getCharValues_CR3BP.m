function [mu,m_star,t_star,l_star,v_star,n_star] = getCharValues_CR3BP(m1,m2,r12)
% =======================================================================
%                      Characteristic Units of CR3BP
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 10-10-2020
%
% Format : [mu,m_star,t_star,l_star,v_star,n_rot] = ...
%               getCharValues_CR3BP(m1,m2,r12)
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% m1            : Mass of primary [1x1]
% m2            : Mass of secondary [1x1]
% r12           : Distance between primary and secondary [1x1]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% mu            : 3-body constant [1x1]
% m_star        : Characteristic Mass [1x1]
% t_star        : Characteristic Time [1x1]
% l_star        : Characteristic Length [1x1]
% v_star        : Characteristic Velocity [1x1]
% n_star        : Characteristic Angular Velocity [1x1]
% -----------------------------------------------------------------------
% 
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 10-10-2020 : Code Created
% -----------------------------------------------------------------------

G = phy_const("C01"); % universal gravitational constant

m_star = m1+m2;

mu = m2/m_star;

l_star = r12;

t_star = sqrt(l_star^3/(G*m_star));

v_star = l_star/t_star;

n_star = sqrt(G*(m1+m2)/r12^3);

end