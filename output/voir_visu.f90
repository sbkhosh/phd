program visu

  implicit none

  integer, parameter :: nx=128,ny=129,nz=128
  integer :: i,j,k,count,stat
  integer :: fn_velx,fn_opac_velx,fn_vely,fn_opac_vely,fn_velz,fn_opac_velz
  integer :: fn_vortx,fn_opac_vortx,fn_vorty,fn_opac_vorty,fn_vortz,fn_opac_vortz
  real(8),dimension(nx) :: y1
  real(8),dimension(ny) :: y2,yp
  real(8),dimension(nz) :: y3
  real(8) :: pi,xlx,zlz,re_tau0,re_p
  character(len=200) :: fiter_velx,fiter_vely,fiter_velz
  character(len=200) :: fiter_vortx,fiter_vorty,fiter_vortz
  real(8), allocatable, dimension(:,:,:) :: velx,vely,velz
  real(8), allocatable, dimension(:,:,:) :: vortx,vorty,vortz

  allocate(velx(nx,ny,nz),vely(nx,ny,nz),velz(nx,ny,nz))
  allocate(vortx(nx,ny,nz),vorty(nx,ny,nz),vortz(nx,ny,nz))

  pi=acos(-1.)
  xlx=4*pi
  zlz=(4.0*pi)/3.0

  re_tau0=179.5
  re_p=4200.0

  open(12,file='ygrid',status='old',iostat=stat)
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

!###############################################################################################################################################
  write(fiter_velx,'(a)') 'velx.dat'
     call get_unit_nr(fn_velx)
     open(fn_velx,file=trim(fiter_velx),form='unformatted',&
          access='direct',recl=8,status='old')
     count = 1
     do k=1,nz
        do j=1,ny
           do i=1,nx
              read(fn_velx,rec=count) velx(i,j,k)
              count = count + 1
           end do
        end do
     end do
     close(fn_velx)
 
  call get_unit_nr(fn_opac_velx)
  open(fn_opac_velx,file='velx.vtr')
  write(fn_opac_velx,*)'<VTKFile type="RectilinearGrid" version="0.1"',' byte_order="LittleEndian">'
  write(fn_opac_velx,*)'  <RectilinearGrid WholeExtent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(fn_opac_velx,*)'    <Piece Extent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(fn_opac_velx,*)'      <Coordinates>'
  write(fn_opac_velx,*)'        <DataArray type="Float64"',' Name="X_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_velx,*) (y1(i),i=1,nx)
  write(fn_opac_velx,*)'        </DataArray>'
  write(fn_opac_velx,*)'        <DataArray type="Float64"',' Name="Y_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_velx,*) (y2(j),j=1,ny)
  write(fn_opac_velx,*)'        </DataArray>'
  write(fn_opac_velx,*)'        <DataArray type="Float64"',' Name="Z_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_velx,*) (y3(k),k=1,nz)
  write(fn_opac_velx,*)'        </DataArray>'
  write(fn_opac_velx,*)'      </Coordinates>'
  write(fn_opac_velx,*)'      <PointData Scalars="scalar">'
  write(fn_opac_velx,*)'        <DataArray Name="test"',' type="Float64"',' NumberOfComponents="1"',' format="ascii">'
  do k=1,nz
     do j=1,ny
        do i=1,nx
           write(fn_opac_velx,*) velx(i,j,k)
        end do
     end do
  end do
  write(fn_opac_velx,*)'        </DataArray>'
  write(fn_opac_velx,*)'      </PointData>'
  write(fn_opac_velx,*)'    </Piece>'
  write(fn_opac_velx,*)'  </RectilinearGrid>'
  write(fn_opac_velx,*)'</VTKFile>'
  close(fn_opac_velx)

!###############################################################################################################################################
  write(fiter_vely,'(a)') 'vely.dat'
     call get_unit_nr(fn_vely)
     open(fn_vely,file=trim(fiter_vely),form='unformatted',&
          access='direct',recl=8,status='old')
     count = 1
     do k=1,nz
        do j=1,ny
           do i=1,nx
              read(fn_vely,rec=count) vely(i,j,k)
              count = count + 1
           end do
        end do
     end do
     close(fn_vely)
 
  call get_unit_nr(fn_opac_vely)
  open(fn_opac_vely,file='vely.vtr')
  write(fn_opac_vely,*)'<VTKFile type="RectilinearGrid" version="0.1"',' byte_order="LittleEndian">'
  write(fn_opac_vely,*)'  <RectilinearGrid WholeExtent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(fn_opac_vely,*)'    <Piece Extent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(fn_opac_vely,*)'      <Coordinates>'
  write(fn_opac_vely,*)'        <DataArray type="Float64"',' Name="X_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_vely,*) (y1(i),i=1,nx)
  write(fn_opac_vely,*)'        </DataArray>'
  write(fn_opac_vely,*)'        <DataArray type="Float64"',' Name="Y_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_vely,*) (y2(j),j=1,ny)
  write(fn_opac_vely,*)'        </DataArray>'
  write(fn_opac_vely,*)'        <DataArray type="Float64"',' Name="Z_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_vely,*) (y3(k),k=1,nz)
  write(fn_opac_vely,*)'        </DataArray>'
  write(fn_opac_vely,*)'      </Coordinates>'
  write(fn_opac_vely,*)'      <PointData Scalars="scalar">'
  write(fn_opac_vely,*)'        <DataArray Name="test"',' type="Float64"',' NumberOfComponents="1"',' format="ascii">'
  do k=1,nz
     do j=1,ny
        do i=1,nx
           write(fn_opac_vely,*) vely(i,j,k)
        end do
     end do
  end do
  write(fn_opac_vely,*)'        </DataArray>'
  write(fn_opac_vely,*)'      </PointData>'
  write(fn_opac_vely,*)'    </Piece>'
  write(fn_opac_vely,*)'  </RectilinearGrid>'
  write(fn_opac_vely,*)'</VTKFile>'
  close(fn_opac_vely)

!###############################################################################################################################################
  write(fiter_velz,'(a)') 'velz.dat'
     call get_unit_nr(fn_velz)
     open(fn_velz,file=trim(fiter_velz),form='unformatted',&
          access='direct',recl=8,status='old')
     count = 1
     do k=1,nz
        do j=1,ny
           do i=1,nx
              read(fn_velz,rec=count) velz(i,j,k)
              count = count + 1
           end do
        end do
     end do
     close(fn_velz)
 
  call get_unit_nr(fn_opac_velz)
  open(fn_opac_velz,file='velz.vtr')
  write(fn_opac_velz,*)'<VTKFile type="RectilinearGrid" version="0.1"',' byte_order="LittleEndian">'
  write(fn_opac_velz,*)'  <RectilinearGrid WholeExtent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(fn_opac_velz,*)'    <Piece Extent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(fn_opac_velz,*)'      <Coordinates>'
  write(fn_opac_velz,*)'        <DataArray type="Float64"',' Name="X_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_velz,*) (y1(i),i=1,nx)
  write(fn_opac_velz,*)'        </DataArray>'
  write(fn_opac_velz,*)'        <DataArray type="Float64"',' Name="Y_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_velz,*) (y2(j),j=1,ny)
  write(fn_opac_velz,*)'        </DataArray>'
  write(fn_opac_velz,*)'        <DataArray type="Float64"',' Name="Z_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_velz,*) (y3(k),k=1,nz)
  write(fn_opac_velz,*)'        </DataArray>'
  write(fn_opac_velz,*)'      </Coordinates>'
  write(fn_opac_velz,*)'      <PointData Scalars="scalar">'
  write(fn_opac_velz,*)'        <DataArray Name="test"',' type="Float64"',' NumberOfComponents="1"',' format="ascii">'
  do k=1,nz
     do j=1,ny
        do i=1,nx
           write(fn_opac_velz,*) velz(i,j,k)
        end do
     end do
  end do
  write(fn_opac_velz,*)'        </DataArray>'
  write(fn_opac_velz,*)'      </PointData>'
  write(fn_opac_velz,*)'    </Piece>'
  write(fn_opac_velz,*)'  </RectilinearGrid>'
  write(fn_opac_velz,*)'</VTKFile>'
  close(fn_opac_velz)

!###############################################################################################################################################
!###############################################################################################################################################
!###############################################################################################################################################

  write(fiter_vortx,'(a)') 'vortx.dat'
     call get_unit_nr(fn_vortx)
     open(fn_vortx,file=trim(fiter_vortx),form='unformatted',&
          access='direct',recl=8,status='old')
     count = 1
     do k=1,nz
        do j=1,ny
           do i=1,nx
              read(fn_vortx,rec=count) vortx(i,j,k)
              count = count + 1
           end do
        end do
     end do
     close(fn_vortx)
 
  call get_unit_nr(fn_opac_vortx)
  open(fn_opac_vortx,file='vortx.vtr')
  write(fn_opac_vortx,*)'<VTKFile type="RectilinearGrid" version="0.1"',' byte_order="LittleEndian">'
  write(fn_opac_vortx,*)'  <RectilinearGrid WholeExtent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(fn_opac_vortx,*)'    <Piece Extent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(fn_opac_vortx,*)'      <Coordinates>'
  write(fn_opac_vortx,*)'        <DataArray type="Float64"',' Name="X_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_vortx,*) (y1(i),i=1,nx)
  write(fn_opac_vortx,*)'        </DataArray>'
  write(fn_opac_vortx,*)'        <DataArray type="Float64"',' Name="Y_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_vortx,*) (y2(j),j=1,ny)
  write(fn_opac_vortx,*)'        </DataArray>'
  write(fn_opac_vortx,*)'        <DataArray type="Float64"',' Name="Z_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_vortx,*) (y3(k),k=1,nz)
  write(fn_opac_vortx,*)'        </DataArray>'
  write(fn_opac_vortx,*)'      </Coordinates>'
  write(fn_opac_vortx,*)'      <PointData Scalars="scalar">'
  write(fn_opac_vortx,*)'        <DataArray Name="test"',' type="Float64"',' NumberOfComponents="1"',' format="ascii">'
  do k=1,nz
     do j=1,ny
        do i=1,nx
           write(fn_opac_vortx,*) vortx(i,j,k)
        end do
     end do
  end do
  write(fn_opac_vortx,*)'        </DataArray>'
  write(fn_opac_vortx,*)'      </PointData>'
  write(fn_opac_vortx,*)'    </Piece>'
  write(fn_opac_vortx,*)'  </RectilinearGrid>'
  write(fn_opac_vortx,*)'</VTKFile>'
  close(fn_opac_vortx)

!###############################################################################################################################################
  write(fiter_vorty,'(a)') 'vorty.dat'
     call get_unit_nr(fn_vorty)
     open(fn_vorty,file=trim(fiter_vorty),form='unformatted',&
          access='direct',recl=8,status='old')
     count = 1
     do k=1,nz
        do j=1,ny
           do i=1,nx
              read(fn_vorty,rec=count) vorty(i,j,k)
              count = count + 1
           end do
        end do
     end do
     close(fn_vorty)
 
  call get_unit_nr(fn_opac_vorty)
  open(fn_opac_vorty,file='vorty.vtr')
  write(fn_opac_vorty,*)'<VTKFile type="RectilinearGrid" version="0.1"',' byte_order="LittleEndian">'
  write(fn_opac_vorty,*)'  <RectilinearGrid WholeExtent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(fn_opac_vorty,*)'    <Piece Extent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(fn_opac_vorty,*)'      <Coordinates>'
  write(fn_opac_vorty,*)'        <DataArray type="Float64"',' Name="X_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_vorty,*) (y1(i),i=1,nx)
  write(fn_opac_vorty,*)'        </DataArray>'
  write(fn_opac_vorty,*)'        <DataArray type="Float64"',' Name="Y_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_vorty,*) (y2(j),j=1,ny)
  write(fn_opac_vorty,*)'        </DataArray>'
  write(fn_opac_vorty,*)'        <DataArray type="Float64"',' Name="Z_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_vorty,*) (y3(k),k=1,nz)
  write(fn_opac_vorty,*)'        </DataArray>'
  write(fn_opac_vorty,*)'      </Coordinates>'
  write(fn_opac_vorty,*)'      <PointData Scalars="scalar">'
  write(fn_opac_vorty,*)'        <DataArray Name="test"',' type="Float64"',' NumberOfComponents="1"',' format="ascii">'
  do k=1,nz
     do j=1,ny
        do i=1,nx
           write(fn_opac_vorty,*) vorty(i,j,k)
        end do
     end do
  end do
  write(fn_opac_vorty,*)'        </DataArray>'
  write(fn_opac_vorty,*)'      </PointData>'
  write(fn_opac_vorty,*)'    </Piece>'
  write(fn_opac_vorty,*)'  </RectilinearGrid>'
  write(fn_opac_vorty,*)'</VTKFile>'
  close(fn_opac_vorty)

!###############################################################################################################################################
  write(fiter_vortz,'(a)') 'vortz.dat'
     call get_unit_nr(fn_vortz)
     open(fn_vortz,file=trim(fiter_vortz),form='unformatted',&
          access='direct',recl=8,status='old')
     count = 1
     do k=1,nz
        do j=1,ny
           do i=1,nx
              read(fn_vortz,rec=count) vortz(i,j,k)
              count = count + 1
           end do
        end do
     end do
     close(fn_vortz)
 
  call get_unit_nr(fn_opac_vortz)
  open(fn_opac_vortz,file='vortz.vtr')
  write(fn_opac_vortz,*)'<VTKFile type="RectilinearGrid" version="0.1"',' byte_order="LittleEndian">'
  write(fn_opac_vortz,*)'  <RectilinearGrid WholeExtent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(fn_opac_vortz,*)'    <Piece Extent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(fn_opac_vortz,*)'      <Coordinates>'
  write(fn_opac_vortz,*)'        <DataArray type="Float64"',' Name="X_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_vortz,*) (y1(i),i=1,nx)
  write(fn_opac_vortz,*)'        </DataArray>'
  write(fn_opac_vortz,*)'        <DataArray type="Float64"',' Name="Y_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_vortz,*) (y2(j),j=1,ny)
  write(fn_opac_vortz,*)'        </DataArray>'
  write(fn_opac_vortz,*)'        <DataArray type="Float64"',' Name="Z_COORDINATES"',' NumberOfComponents="1">'
  write(fn_opac_vortz,*) (y3(k),k=1,nz)
  write(fn_opac_vortz,*)'        </DataArray>'
  write(fn_opac_vortz,*)'      </Coordinates>'
  write(fn_opac_vortz,*)'      <PointData Scalars="scalar">'
  write(fn_opac_vortz,*)'        <DataArray Name="test"',' type="Float64"',' NumberOfComponents="1"',' format="ascii">'
  do k=1,nz
     do j=1,ny
        do i=1,nx
           write(fn_opac_vortz,*) vortz(i,j,k)
        end do
     end do
  end do
  write(fn_opac_vortz,*)'        </DataArray>'
  write(fn_opac_vortz,*)'      </PointData>'
  write(fn_opac_vortz,*)'    </Piece>'
  write(fn_opac_vortz,*)'  </RectilinearGrid>'
  write(fn_opac_vortz,*)'</VTKFile>'
  close(fn_opac_vortz)

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
