function Kmstab=rigistab(BXIdev,BETAdev,BXIETAdev,BETAZETAdev,BXIZETAdev,C)
 
 Kmstab=8*(BXIdev'*C*BXIdev+BETAdev'*C*BETAdev)/3+...
        8*(BETAZETAdev'*C*BETAZETAdev+...
           BXIETAdev'*C*BXIETAdev+...
           BXIZETAdev'*C*BXIZETAdev)/9.0;

end