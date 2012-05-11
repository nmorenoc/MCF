  vi_new(1) = &
       ((eij(1)*eij(2)*(vj_old(2)-vi_old(2))+&
       2*vi_old(1)*eij(2)**2+eij(1)**2*vj_old(1)+eij(1)**2*vi_old(1))*Aij-&
       vi_old(1))/((2*eij(2)**2+2*eij(1)**2)*Aij-1)
  vi_new(2) = &
       ((eij(2)**2*(vj_old(2)+vi_old(2))+2*eij(1)**2*vi_old(2)+&
       (eij(1)*vj_old(1)-eij(1)*vi_old(1))*eij(2))*Aij-vi_old(2))/&
       ((2*eij(2)**2+2*eij(1)**2)*Aij-1)
  vj_new(1) = &
       ((eij(1)*eij(2)*(vi_old(2)-vj_old(2))+&
       2*vj_old(1)*eij(2)**2+eij(1)**2*vj_old(1)+eij(1)**2*vi_old(1))*Aij-&
       vj_old(1))/((2*eij(2)**2+2*eij(1)**2)*Aij-1)
  vj_new(2) = &
       ((eij(2)**2*(vj_old(2)+vi_old(2))+2*eij(1)**2*vj_old(2)+&
       (eij(1)*vi_old(1)-eij(1)*vj_old(1))*eij(2))*Aij-vj_old(2))/&
       ((2*eij(2)**2+2*eij(1)**2)*Aij-1)

