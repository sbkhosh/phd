!    int_uv_lw=0.0;int_uv_uw=0.0
!    do j=2,ny
!       int_uv_lw=int_uv_lw+0.5*(((2.0-yp(j+1))*(S_uv(j+1)))+((2.0-yp(j))*(S_uv(j))))*&
!            (yp(j+1)-yp(j))
!       int_uv_lw=int_uv_lw+0.5*(((2.0-yp(j+1))*(S_ubar(j+1)*S_vbar(j+1)))+((2.0-yp(j))*(S_ubar(j)*S_vbar(j))))*&
!            (yp(j+1)-yp(j))

!       int_uv_uw=int_uv_uw+0.5*(((2.0-yp(j+1))*(S_uv(j+1)))+((2.0-yp(j))*(S_uv(j))))*&
!            (yp(j+1)-yp(j))
!       int_uv_uw=int_uv_uw+0.5*(((2.0-yp(j+1))*(S_ubar(j+1)*S_vbar(j+1)))+((2.0-yp(j))*(S_ubar(j)*S_vbar(j))))*&
!            (yp(j+1)-yp(j))
!    end do
!    int_uv_lw=-12.0*int_uv_lw
!    int_uv_lw=int_uv_lw+24.0*(S_ubar(1)*S_vbar(1)+S_uv(1))
!    int_uv_lw=int_uv_lw+16.0*(S_ubar((ny-1)/2)*S_vbar((ny-1)/2)-S_ubar(1)*S_vbar(1)+S_uv((ny-1)/2)-S_uv(1))
!    int_uv_lw=int_uv_lw/(16.0/9.0)
!    int_uv_lw=int_uv_lw+(12.0/(5600.0*3.0/4.0))*(1.0-2.0*S_ubar(1)/(4.0/3.0))

!    int_uv_uw=-12.0*int_uv_uw
!    int_uv_uw=int_uv_uw+24.0*(S_ubar(ny)*S_vbar(ny)+S_uv(ny))
!    int_uv_uw=int_uv_uw+16.0*(S_ubar(ny)*S_vbar(ny)-S_ubar((ny-1)/2)*S_vbar((ny-1)/2)+S_uv(ny)-S_uv((ny-1)/2))
!    int_uv_uw=int_uv_uw/(16.0/9.0)
!    int_uv_uw=int_uv_uw+(12.0/(5600.0*3.0/4.0))*(1.0-2.0*S_ubar(ny)/(4.0/3.0))

!    write(6,*) "#####################################################################"
!    write(6,*) "cf_fik_lw = ", int_uv_lw
!    write(6,*) "cf_fik_uw = ", int_uv_uw
!    write(6,*) "cf_fik_tot = ", (int_uv_lw+int_uv_uw)/2.0
!    write(6,*) "cf_tot_code = ", c_f_tot
!    write(6,*) "#####################################################################"


!    call WRITE_STATISTIC(S_ubar,S_urms,S_vbar,S_vrms,S_wbar,S_wrms,&
!         S_uv,S_uw,S_vw,S_dudy_lw_tot,S_dudy_uw_tot,S_dudy_all_tot,&
!         c_f_tot,S_vortxbar,S_vortxrms,S_vortybar,S_vortyrms,S_vortzbar,S_vortzrms)
!    close(10)



!############################################################################
!
subroutine STATISTIC_WO_VORT(ux1,uy1,uz1,dudy1,ta1,tb1,tc1,td1,te1,tf1,tg1,&
     th1,ti1,di1,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io
USE MPI
USE var, only: S_dudy_lw_tot,S_dudy_uw_tot,S_dudy_tot,c_f_tot,&
     S_ubar,S_vbar,S_wbar,S_urms,S_vrms,S_wrms,S_uv,S_uw,S_vw,&
     S_vortxbar,S_vortybar,S_vortzbar,S_vortxrms,S_vortyrms,S_vortzrms

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,dudy1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tux1,tuy1,tuz1,tdudy1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: vortx1,vorty1,vortz1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tvortx1,tvorty1,tvortz1
real(mytype),dimension(ysize(2)) :: ubar,ubar_red,urms,urms_red
real(mytype),dimension(ysize(2)) :: vbar,vbar_red,vrms,vrms_red
real(mytype),dimension(ysize(2)) :: wbar,wbar_red,wrms,wrms_red
real(mytype),dimension(ysize(2)) :: vortxbar,vortxbar_red,vortxrms,vortxrms_red
real(mytype),dimension(ysize(2)) :: vortybar,vortybar_red,vortyrms,vortyrms_red
real(mytype),dimension(ysize(2)) :: vortzbar,vortzbar_red,vortzrms,vortzrms_red
real(mytype),dimension(ysize(2)) :: uv,uv_red,uw,uw_red,vw,vw_red
real(mytype) :: dudy_mean_lw,dudy_mean_uw,dudy_mean,c_f
real(mytype) :: dudy_mean_lw_red,dudy_mean_uw_red
real(mytype),dimension(5,ilast-ifirst+1) :: drg_inst
integer :: i,j,k,ss,code
integer, dimension(2) :: dims, dummy_coords
logical, dimension(2) :: dummy_periods

TYPE(DECOMP_INFO) :: phG
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3

if(itime.eq.ifirst.and.nrank.eq.0) then
   open(500,file='drg.txt',access='append')
end if

nvect1=xsize(1)*xsize(2)*xsize(3)
!x-derivatives
call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!y-derivatives
call transpose_x_to_y(ux1,td2)
call transpose_x_to_y(uy1,te2)
call transpose_x_to_y(uz1,tf2)
call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!!z-derivatives
call transpose_y_to_z(td2,td3)
call transpose_y_to_z(te2,te3)
call transpose_y_to_z(tf2,tf3)
call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
!!all back to x-pencils
call transpose_z_to_y(ta3,td2)
call transpose_z_to_y(tb3,te2)
call transpose_z_to_y(tc3,tf2)
call transpose_y_to_x(td2,tg1)
call transpose_y_to_x(te2,th1)
call transpose_y_to_x(tf2,ti1)
call transpose_y_to_x(ta2,td1)
call transpose_y_to_x(tb2,te1)
call transpose_y_to_x(tc2,tf1)

!du/dx=ta1 du/dy=td1 and du/dz=tg1
!dv/dx=tb1 dv/dy=te1 and dv/dz=th1
!dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

!############################################################################
!VORTICITY - x component
vortxtf1(ijk,1,1)-th1(ijk,1,1)
!############################################################################
!VORTICITY - y component
tg1(ijk,1,1)-tc1(ijk,1,1)
!############################################################################
!VORTICITY - z component
di1(ijk,1,1)=tb1(ijk,1,1)-td1(ijk,1,1)
!############################################################################

call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
     dims, dummy_periods, dummy_coords, code)

!!!!!!!!!!!!!!!!!!!
! mean quantities !
!!!!!!!!!!!!!!!!!!!
do j=1,ny
   ubar(j)=0.0;vbar(j)=0.0;wbar(j)=0.0
   vortxbar(j)=0.0;vortybar(j)=0.0;vortzbar(j)=0.0
end do

call transpose_x_to_y(ux1,tux1)
call transpose_x_to_y(uy1,tuy1)
call transpose_x_to_y(uz1,tuz1)

call transpose_x_to_y(vortx1,tvortx1)
call transpose_x_to_y(vorty1,tvorty1)
call transpose_x_to_y(vortz1,tvortz1)

do j=1,ny
   do k=1,ysize(3)
      do i=1,ysize(1)
         ubar(j)=ubar(j)+tux1(i,j,k)
         vbar(j)=vbar(j)+tuy1(i,j,k)
         wbar(j)=wbar(j)+tuz1(i,j,k)
         vortxbar(j)=vortxbar(j)+tvortx1(i,j,k)
         vortybar(j)=vortybar(j)+tvorty1(i,j,k)
         vortzbar(j)=vortzbar(j)+tvortz1(i,j,k)
      end do
   end do
end do

do j=1,ny
   ubar(j)=ubar(j)/ysize(1)/ysize(3)
   vbar(j)=vbar(j)/ysize(1)/ysize(3)
   wbar(j)=wbar(j)/ysize(1)/ysize(3)
   vortxbar(j)=vortxbar(j)/ysize(1)/ysize(3)
   vortybar(j)=vortybar(j)/ysize(1)/ysize(3)
   vortzbar(j)=vortzbar(j)/ysize(1)/ysize(3)
end do

call MPI_ALLREDUCE(ubar,ubar_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(vbar,vbar_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(wbar,wbar_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(vortxbar,vortxbar_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(vortybar,vortybar_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(vortzbar,vortzbar_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)

do j=1,ny
   ubar_red(j)=ubar_red(j)/nproc
   vbar_red(j)=vbar_red(j)/nproc
   wbar_red(j)=wbar_red(j)/nproc
   vortxbar_red(j)=vortxbar_red(j)/nproc
   vortybar_red(j)=vortybar_red(j)/nproc
   vortzbar_red(j)=vortzbar_red(j)/nproc
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! rms quantities and Reynolds stresses !            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do j=1,ny
   urms(j)=0.0;vrms(j)=0.0;wrms(j)=0.0
   uv(j)=0.0;uw(j)=0.0;vw(j)=0.0
   vortxrms(j)=0.0;vortyrms(j)=0.0;vortzrms(j)=0.0
end do

do j=1,ny
   do k=1,ysize(3)
      do i=1,ysize(1)
         urms(j)=urms(j)+(tux1(i,j,k)-ubar_red(j))**2
         vrms(j)=vrms(j)+(tuy1(i,j,k)-vbar_red(j))**2
         wrms(j)=wrms(j)+(tuz1(i,j,k)-wbar_red(j))**2
         uv(j)=uv(j)+(tux1(i,j,k)-ubar_red(j))*(tuy1(i,j,k)-vbar_red(j))
         uw(j)=uw(j)+(tux1(i,j,k)-ubar_red(j))*(tuz1(i,j,k)-wbar_red(j))
         vw(j)=vw(j)+(tuy1(i,j,k)-vbar_red(j))*(tuz1(i,j,k)-wbar_red(j))
         vortxrms(j)=vortxrms(j)+(tvortx1(i,j,k)-vortxbar_red(j))**2
         vortyrms(j)=vortyrms(j)+(tvorty1(i,j,k)-vortybar_red(j))**2
         vortzrms(j)=vortzrms(j)+(tvortz1(i,j,k)-vortzbar_red(j))**2
      end do
   end do
end do

do j=1,ny
   urms(j)=urms(j)/ysize(1)/ysize(3)
   vrms(j)=vrms(j)/ysize(1)/ysize(3)
   wrms(j)=wrms(j)/ysize(1)/ysize(3)
   uv(j)=uv(j)/ysize(1)/ysize(3)
   uw(j)=uw(j)/ysize(1)/ysize(3)
   vw(j)=vw(j)/ysize(1)/ysize(3)
   vortxrms(j)=vortxrms(j)/ysize(1)/ysize(3)
   vortyrms(j)=vortyrms(j)/ysize(1)/ysize(3)
   vortzrms(j)=vortzrms(j)/ysize(1)/ysize(3)
end do

call MPI_ALLREDUCE(urms,urms_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(vrms,vrms_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(wrms,wrms_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(uv,uv_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(uw,uw_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(vw,vw_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(vortxrms,vortxrms_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(vortyrms,vortyrms_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(vortzrms,vortzrms_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)

do j=1,ny
   urms_red(j)=urms_red(j)/nproc
   vrms_red(j)=vrms_red(j)/nproc
   wrms_red(j)=wrms_red(j)/nproc
   uv_red(j)=uv_red(j)/nproc
   uw_red(j)=uw_red(j)/nproc
   vw_red(j)=vw_red(j)/nproc
   vortxrms_red(j)=vortxrms_red(j)/nproc
   vortyrms_red(j)=vortyrms_red(j)/nproc
   vortzrms_red(j)=vortzrms_red(j)/nproc
end do

do j=1,ny
   S_ubar(j)=S_ubar(j)+ubar_red(j)
   S_vbar(j)=S_vbar(j)+vbar_red(j)
   S_wbar(j)=S_wbar(j)+wbar_red(j)
   S_urms(j)=S_urms(j)+urms_red(j)
   S_vrms(j)=S_vrms(j)+vrms_red(j)
   S_wrms(j)=S_wrms(j)+wrms_red(j)
   S_uv(j)=S_uv(j)+uv_red(j)
   S_uw(j)=S_uw(j)+uw_red(j)
   S_vw(j)=S_vw(j)+vw_red(j)
   S_vortxrms(j)=S_vortxrms(j)+vortxrms_red(j)
   S_vortyrms(j)=S_vortyrms(j)+vortyrms_red(j)
   S_vortzrms(j)=S_vortzrms(j)+vortzrms_red(j)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! drag-related quantities !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

call transpose_x_to_y(dudy1,tdudy1)
dudy_mean_lw=0.0
dudy_mean_uw=0.0

do k=1,ysize(3)
   do i=1,ysize(1)
      dudy_mean_lw=dudy_mean_lw+tdudy1(i,1,k)
      dudy_mean_uw=dudy_mean_uw+tdudy1(i,ny,k)
   end do
end do
dudy_mean_lw=dudy_mean_lw/ysize(1)/ysize(3)
dudy_mean_uw=dudy_mean_uw/ysize(1)/ysize(3)

call mpi_allreduce(dudy_mean_lw,dudy_mean_lw_red,1,real_type,&
     MPI_SUM,MPI_COMM_WORLD,code)
call mpi_allreduce(dudy_mean_uw,dudy_mean_uw_red,1,real_type,&
     MPI_SUM,MPI_COMM_WORLD,code)

dudy_mean_lw_red=dudy_mean_lw_red/nproc
dudy_mean_uw_red=dudy_mean_uw_red/nproc

dudy_mean=(dudy_mean_lw_red-dudy_mean_uw_red)/(2.0)
c_f=4.5*dudy_mean/re

if(nrank.eq.0) then
   ss=itime+(ilast-ifirst+1)-ilast
   drg_inst(1,ss)=dble(itime*imodulo*dt)
   drg_inst(2,ss)=dudy_mean_lw_red
   drg_inst(3,ss)=dudy_mean_uw_red
   drg_inst(4,ss)=dudy_mean
   drg_inst(5,ss)=c_f
   write(500,*) drg_inst(1,ss),drg_inst(2,ss),drg_inst(3,ss),&
        drg_inst(4,ss),drg_inst(5,ss)
end if

S_dudy_lw_tot=S_dudy_lw_tot+dudy_mean_lw
S_dudy_uw_tot=S_dudy_uw_tot+dudy_mean_uw
S_dudy_tot=S_dudy_tot+dudy_mean
c_f_tot=c_f_tot+c_f

if(itime.eq.ilast.and.nrank.eq.0) then
   do j=1,ny
      S_ubar(j)=S_ubar(j)/(dble(ilast-ifirst+1))
      S_vbar(j)=S_vbar(j)/(dble(ilast-ifirst+1))
      S_wbar(j)=S_wbar(j)/(dble(ilast-ifirst+1))
      S_urms(j)=S_urms(j)/(dble(ilast-ifirst+1))
      S_urms(j)=sqrt(S_urms(j))
      S_vrms(j)=S_vrms(j)/(dble(ilast-ifirst+1))
      S_vrms(j)=sqrt(S_vrms(j))
      S_wrms(j)=S_wrms(j)/(dble(ilast-ifirst+1))
      S_wrms(j)=sqrt(S_wrms(j))
      S_uv(j)=S_uv(j)/(dble(ilast-ifirst+1))
      S_uw(j)=S_uw(j)/(dble(ilast-ifirst+1))
      S_vw(j)=S_vw(j)/(dble(ilast-ifirst+1))
      S_vortxrms(j)=S_vortxrms(j)/(dble(ilast-ifirst+1))
      S_vortxrms(j)=sqrt(S_vortxrms(j))
      S_vortyrms(j)=S_vortyrms(j)/(dble(ilast-ifirst+1))
      S_vortyrms(j)=sqrt(S_vortyrms(j))
      S_vortzrms(j)=S_vortzrms(j)/(dble(ilast-ifirst+1))
      S_vortzrms(j)=sqrt(S_vortzrms(j))
   end do

   S_dudy_lw_tot=S_dudy_lw_tot/(dble(ilast-ifirst+1))
   S_dudy_uw_tot=S_dudy_uw_tot/(dble(ilast-ifirst+1))
   S_dudy_tot=S_dudy_tot/(dble(ilast-ifirst+1))
   c_f_tot=c_f_tot/(dble(ilast-ifirst+1))
end if

end subroutine STATISTIC_WO_VORT


!############################################################################
!
subroutine CALC_VORT(ux1,uy1,uz1,ta1,tb1,tc1,td1,te1,tf1,tg1,&
     th1,ti1,di1,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,vrtx,vrty,vrtz)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io
USE MPI

implicit none

TYPE(DECOMP_INFO) :: phG
integer :: i,j,k
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(out) :: vrtx,vrty,vrtz

!x-derivatives
call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!y-derivatives
call transpose_x_to_y(ux1,td2)
call transpose_x_to_y(uy1,te2)
call transpose_x_to_y(uz1,tf2)
call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!!z-derivatives
call transpose_y_to_z(td2,td3)
call transpose_y_to_z(te2,te3)
call transpose_y_to_z(tf2,tf3)
call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
!!all back to x-pencils
call transpose_z_to_y(ta3,td2)
call transpose_z_to_y(tb3,te2)
call transpose_z_to_y(tc3,tf2)
call transpose_y_to_x(td2,tg1)
call transpose_y_to_x(te2,th1)
call transpose_y_to_x(tf2,ti1)
call transpose_y_to_x(ta2,td1)
call transpose_y_to_x(tb2,te1)
call transpose_y_to_x(tc2,tf1)

!du/dx=ta1 du/dy=td1 and du/dz=tg1
!dv/dx=tb1 dv/dy=te1 and dv/dz=th1
!dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

!############################################################################
!VORTICITY - x-y-z component
do k=1,xsize(3)
   do j=1,xsize(2)
      do i=1,xsize(1)
         vrtx(i,j,k)=tf1(i,j,k)-th1(i,j,k)
         vrty(i,j,k)=tg1(i,j,k)-tc1(i,j,k)
         vrtz(i,j,k)=tb1(i,j,k)-td1(i,j,k)
      end do
   end do
end do

return
end subroutine CALC_VORT



     if(with_vort_dat) then
        write(filename, "(a,i6.6,a)") trim(datdir)//"/vortx", itime ,".dat"
        call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(filename), &
             MPI_MODE_RDONLY, MPI_INFO_NULL, &
             fhvox, ierrorvox)
        call decomp_2d_read_var(fhvox,dispvox,1,tabvox)
        call MPI_FILE_CLOSE(fhvox,ierrorvox)
        
        write(filename, "(a,i6.6,a)") trim(datdir)//"/vorty", itime ,".dat"
        call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(filename), &
             MPI_MODE_RDONLY, MPI_INFO_NULL, &
             fhvoy, ierrorvoy)
        call decomp_2d_read_var(fhvoy,dispvoy,1,tabvoy)
        call MPI_FILE_CLOSE(fhvoy,ierrorvoy)
        
        write(filename, "(a,i6.6,a)") trim(datdir)//"/vortz", itime ,".dat"
        call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(filename), &
             MPI_MODE_RDONLY, MPI_INFO_NULL, &
             fhvoz, ierrorvoz)
        call decomp_2d_read_var(fhvoz,dispvoz,1,tabvoz)
        call MPI_FILE_CLOSE(fhvoz,ierrorvoz)
     end if

     ! write(filename, "(a,i6.6,a)") trim(datdir)//"/dudy", itime ,".dat"
     ! call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(filename), &
     !      MPI_MODE_RDONLY, MPI_INFO_NULL, &
     !      fhdudy, ierrordudy)
     ! call decomp_2d_read_var(fhdudy,dispdudy,1,tabdudy)
     ! call MPI_FILE_CLOSE(fhdudy,ierrordudy)
     



!############################################################################
!
subroutine CONV_TRSP(ux1,uy1,uz1)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io
USE FR_perf_mod
USE MPI
USE var, only: conv_trspx123_tavg,conv_trspy123_tavg,conv_trspz123_tavg

implicit none

integer :: i,j,k
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: conv_trspx1,conv_trspx2,conv_trspx3,di1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: dconv_trspx1,dconv_trspx2,dconv_trspx3
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tdconv_trspx123

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: conv_trspy1,conv_trspy2,conv_trspy3,di2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tconv_trspy1,tconv_trspy2,tconv_trspy3
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: dtconv_trspy1,dtconv_trspy2,dtconv_trspy3

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: conv_trspz1,conv_trspz2,conv_trspz3,di3
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tconv_trspz1,tconv_trspz2,tconv_trspz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ttconv_trspz1,ttconv_trspz2,ttconv_trspz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: dttconv_trspz1,dttconv_trspz2,dttconv_trspz3
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tdconv_trspz123

real(mytype),dimension(ysize(2)) :: spdconv_trspx123,spdconv_trspy123,spdconv_trspz123

call perfon('conv_trsp')

do k=1,xsize(3)
   do j=1,xsize(2)
      do i=1,xsize(1)
         conv_trspx1(i,j,k)=ux1(i,j,k)*(0.5*ux1(i,j,k)*ux1(i,j,k))
         conv_trspx2(i,j,k)=ux1(i,j,k)*(0.5*uy1(i,j,k)*uy1(i,j,k))
         conv_trspx3(i,j,k)=ux1(i,j,k)*(0.5*uz1(i,j,k)*uz1(i,j,k))

         conv_trspy1(i,j,k)=uy1(i,j,k)*(0.5*ux1(i,j,k)*ux1(i,j,k))
         conv_trspy2(i,j,k)=uy1(i,j,k)*(0.5*uy1(i,j,k)*uy1(i,j,k))
         conv_trspy3(i,j,k)=uy1(i,j,k)*(0.5*uz1(i,j,k)*uz1(i,j,k))

         conv_trspz1(i,j,k)=uz1(i,j,k)*(0.5*ux1(i,j,k)*ux1(i,j,k))
         conv_trspz2(i,j,k)=uz1(i,j,k)*(0.5*uy1(i,j,k)*uy1(i,j,k))
         conv_trspz3(i,j,k)=uz1(i,j,k)*(0.5*uz1(i,j,k)*uz1(i,j,k))
      end do
   end do
end do

call derx (dconv_trspx1,conv_trspx1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call derx (dconv_trspx2,conv_trspx2,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call derx (dconv_trspx3,conv_trspx3,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)

call transpose_x_to_y(conv_trspy1,tconv_trspy1)
call transpose_x_to_y(conv_trspy2,tconv_trspy2)
call transpose_x_to_y(conv_trspy2,tconv_trspy3)
call dery (dtconv_trspy1,tconv_trspy1,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call dery (dtconv_trspy2,tconv_trspy2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call dery (dtconv_trspy3,tconv_trspy3,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)

call transpose_x_to_y(conv_trspz1,tconv_trspz1)
call transpose_x_to_y(conv_trspz2,tconv_trspz2)
call transpose_x_to_y(conv_trspz2,tconv_trspz3)
call transpose_y_to_z(tconv_trspz1,ttconv_trspz1)
call transpose_y_to_z(tconv_trspz2,ttconv_trspz2)
call transpose_y_to_z(tconv_trspz2,ttconv_trspz3)
call derz (dttconv_trspz1,ttconv_trspz1,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (dttconv_trspz2,ttconv_trspz2,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (dttconv_trspz3,ttconv_trspz3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)

call transpose_x_to_y(dconv_trspx1+dconv_trspx2+dconv_trspx3,tdconv_trspx123)
call space_avg(tdconv_trspx123,spdconv_trspx123)

call space_avg(dtconv_trspy1+dtconv_trspy2+dtconv_trspy3,spdconv_trspy123)

call transpose_z_to_y(dttconv_trspz1+dttconv_trspz2+dttconv_trspz3,tdconv_trspz123)
call space_avg(tdconv_trspz123,spdconv_trspz123)

do j=1,ny
   conv_trspx123_tavg(j)=conv_trspx123_tavg(j)+spdconv_trspx123(j)
   conv_trspy123_tavg(j)=conv_trspy123_tavg(j)+spdconv_trspy123(j)
   conv_trspz123_tavg(j)=conv_trspz123_tavg(j)+spdconv_trspz123(j)
end do

call perfoff
end subroutine CONV_TRSP

         ! (2.0*ddudx1(i,j,k))*ddudx1(i,j,k)+&
         !      (ddudy1(i,j,k)+ddvdx1(i,j,k))*ddudy1(i,j,k)+&
         !      (ddudz1(i,j,k)+ddwdx1(i,j,k))*ddudz1(i,j,k)+&
         !      (ddvdx1(i,j,k)+ddudy1(i,j,k))*ddvdx1(i,j,k)+&
         !      (2.0*ddvdy1(i,j,k))*ddvdy1(i,j,k)+&
         !      (ddvdz1(i,j,k)+ddwdy1(i,j,k))*ddvdz1(i,j,k)+&
         !      (ddwdx1(i,j,k)+ddudz1(i,j,k))*ddwdx1(i,j,k)+&
         !      (ddwdy1(i,j,k)+ddvdz1(i,j,k))*ddwdy1(i,j,k)+&
         !      (2.0*ddwdz1(i,j,k))*ddwdz1(i,j,k)


!############################################################################
!
subroutine NRG_SPECTRA(vel,fft_vel)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io
USE FR_perf_mod
USE MPI
USE var, only: nrg_uu_fft,nrg_vv_fft,nrg_ww_fft

implicit none

integer :: i,j,k
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,uuv,uv,vvv,wwv

call perfon('nrg_spectra')

call perfoff
end subroutine NRG_SPECTRA
