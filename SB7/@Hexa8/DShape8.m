function [dN_dXsi,dN_dEta,dN_dZeta]=DShape8(obj,Xsi,Eta,Zeta)
% dérivées des fonctions de formes par rapport au repère isoparamétrique

dN_dXsi(1)=0.125*(-1)*(1-Eta)*(1-Zeta); % noeud 1
dN_dXsi(2)=0.125*(+1)*(1-Eta)*(1-Zeta); % noeud 2 
dN_dXsi(3)=0.125*(+1)*(1+Eta)*(1-Zeta); % noeud 3 
dN_dXsi(4)=0.125*(-1)*(1+Eta)*(1-Zeta); % noeud 4 
dN_dXsi(5)=0.125*(-1)*(1-Eta)*(1+Zeta); % noeud 5 
dN_dXsi(6)=0.125*(+1)*(1-Eta)*(1+Zeta); % noeud 6 
dN_dXsi(7)=0.125*(+1)*(1+Eta)*(1+Zeta); % noeud 7 
dN_dXsi(8)=0.125*(-1)*(1+Eta)*(1+Zeta); % noeud 8 


dN_dEta(1)=0.125*(1-Xsi)*(-1)*(1-Zeta); % noeud 1
dN_dEta(2)=0.125*(1+Xsi)*(-1)*(1-Zeta); % noeud 2
dN_dEta(3)=0.125*(1+Xsi)*(+1)*(1-Zeta); % noeud 3
dN_dEta(4)=0.125*(1-Xsi)*(+1)*(1-Zeta); % noeud 4
dN_dEta(5)=0.125*(1-Xsi)*(-1)*(1+Zeta); % noeud 5
dN_dEta(6)=0.125*(1+Xsi)*(-1)*(1+Zeta); % noeud 6
dN_dEta(7)=0.125*(1+Xsi)*(+1)*(1+Zeta); % noeud 7
dN_dEta(8)=0.125*(1-Xsi)*(+1)*(1+Zeta); % noeud 8


dN_dZeta(1)=0.125*(1-Xsi)*(1-Eta)*(-1); % noeud 1
dN_dZeta(2)=0.125*(1+Xsi)*(1-Eta)*(-1); % noeud 2
dN_dZeta(3)=0.125*(1+Xsi)*(1+Eta)*(-1); % noeud 3
dN_dZeta(4)=0.125*(1-Xsi)*(1+Eta)*(-1); % noeud 4
dN_dZeta(5)=0.125*(1-Xsi)*(1-Eta)*(+1); % noeud 5
dN_dZeta(6)=0.125*(1+Xsi)*(1-Eta)*(+1); % noeud 6
dN_dZeta(7)=0.125*(1+Xsi)*(1+Eta)*(+1); % noeud 7
dN_dZeta(8)=0.125*(1-Xsi)*(1+Eta)*(+1); % noeud 8
end