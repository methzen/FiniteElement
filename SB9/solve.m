function S=solve(GKF,RESIDU,SDISPT,NEQ,ITER,PROP)
      % Prescribed displacement BC
%             size(GKF)
%       size(RESIDU)

alpha=PROP(1);
      NDISP=size(SDISPT,1);
      SDISP = SDISPT(:,3);
      if NDISP~=0
       FIXEDDOF=3*(SDISPT(:,1)-1)+SDISPT(:,2); % global connect
       GKF(FIXEDDOF,:)=zeros(NDISP,NEQ);
       GKF(FIXEDDOF,FIXEDDOF)=alpha*eye(NDISP);
       %
       RESIDU(FIXEDDOF)=0;
       if ITER==1, RESIDU(FIXEDDOF) =alpha*SDISP(:); end
      end
      S=GKF\RESIDU;
end