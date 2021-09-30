function RESN=checkconvergence(FORCE,SDISPT,NEQ,FORCEXTG)    
% Check convergence
        FIXEDDOF=3*(SDISPT(:,1)-1)+SDISPT(:,2);
        ALLDOF=1:NEQ;
        FREEDOF=setdiff(ALLDOF,FIXEDDOF); % take the free dof
        RESN=max(abs(FORCE(FREEDOF)))/max(abs(FORCEXTG));
        RESN=full(RESN);
        %
end