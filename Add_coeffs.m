function  [ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,ind_add,val_add)
%  Add  coefficients  to  an  array
%
if  (length(val_add(:))==1)
val_add   =     ones(size(ind_add))*val_add;
end
ind_list       =     [ind_list,  ind_add(:)];
ind_val        =     [ind_val ,  val_add(:)];