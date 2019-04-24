program visu

  implicit none

  integer, parameter :: nx=256,ny=129,nz=256
  real(8),dimension(nx,ny,nz) :: uxd_bin_file
  real(8),dimension(nx,ny,nz) :: uzd_bin_file
  real(8),dimension(nx,ny,nz) :: uplane_bin_file
  integer :: i,j,k,count,nfilx,nfilz,nfilsqr,nfilplane,stat
  integer :: fnx,fnz,fnplane
  real(8),dimension(nx) :: y1
  real(8),dimension(ny) :: y2,yp
  real(8),dimension(nz) :: y3
  real(8) :: pi,xlx,zlz
  character(len=128) :: fiter_x,fiter_z,fiter_plane
  character(len=128) :: fiter_tampon_x,fiter_tampon_z,fiter_tampon_sqr,fiter_tampon_plane

  pi=acos(-1.)
  xlx=4.53*pi
  zlz=2.26*pi

  open(12,file='yp.dat',status='old',iostat=stat)
  read(12,'(f20.12)') (yp(j), j=1,ny)
  close(12)

  do i=1,nx
     y1(i)=(i-1)*real(xlx)/real(nx)
  enddo

  do j=1,ny
     y2(j)=yp(j)
  enddo

  do k=1,nz
     y3(k)=(k-1)*real(zlz)/real(nz)
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! writing constant plane velocity in vtk format !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(fiter_plane,'(a)') 'disks_opac_plane.dat'

  write(fiter_tampon_plane,'(a)') 'tampon_opac_plane.vtr'
  call get_unit_nr(fnplane)
  open(fnplane,file=trim(fiter_plane),form='unformatted',&
       access='direct',recl=8,status='old')
  count = 1
  do k=1,nz
     do j=1,ny
        do i=1,nx
           read(fnplane,rec=count) uplane_bin_file(i,j,k)
           count = count + 1
        end do
     end do
  end do
  close(fnplane)

  call get_unit_nr(nfilplane)
  open(nfilplane,file=trim(fiter_tampon_plane))
  write(nfilplane,*)'<VTKFile type="RectilinearGrid" version="0.1"',' byte_order="LittleEndian">'
  write(nfilplane,*)'  <RectilinearGrid WholeExtent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(nfilplane,*)'    <Piece Extent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(nfilplane,*)'      <Coordinates>'
  write(nfilplane,*)'        <DataArray type="Float64"',' Name="X_COORDINATES"',' NumberOfComponents="1">'
  write(nfilplane,*) (y1(i),i=1,nx)
  write(nfilplane,*)'        </DataArray>'
  write(nfilplane,*)'        <DataArray type="Float64"',' Name="Y_COORDINATES"',' NumberOfComponents="1">'
  write(nfilplane,*) (y2(j),j=1,ny)
  write(nfilplane,*)'        </DataArray>'
  write(nfilplane,*)'        <DataArray type="Float64"',' Name="Z_COORDINATES"',' NumberOfComponents="1">'
  write(nfilplane,*) (y3(k),k=1,nz)
  write(nfilplane,*)'        </DataArray>'
  write(nfilplane,*)'      </Coordinates>'
  write(nfilplane,*)'      <PointData Scalars="scalar">'
  write(nfilplane,*)'        <DataArray Name="test"',' type="Float64"',' NumberOfComponents="1"',' format="ascii">'
  write(nfilplane,*) (((uplane_bin_file(i,j,k),i=1,nx),j=1,ny),k=1,nz)
  write(nfilplane,*)'        </DataArray>'
  write(nfilplane,*)'      </PointData>'
  write(nfilplane,*)'    </Piece>'
  write(nfilplane,*)'  </RectilinearGrid>'
  write(nfilplane,*)'</VTKFile>'
  close(nfilplane)
end program visu

subroutine get_unit_nr(number)

  integer,intent(out):: number
  ! local variables
  logical :: already_used
  number=31
  do
     inquire(unit=number,opened=already_used)
     if (already_used) then
        number = number + 1
     else
        ! found a free file unit
        exit
     end if
  end do
end subroutine get_unit_nr
