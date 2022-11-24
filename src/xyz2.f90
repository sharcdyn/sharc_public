program xyz2geom
implicit none
integer nat
character*2, allocatable :: at(:)
real*8, allocatable :: geom(:,:),vel(:,:)
real*8, allocatable :: noat(:),mass(:),chat(:)
character*1000 kkchar
character*500 kk

real*8 :: bohr2ang=0.52917721067

integer i,j,k,n,l
integer nfich

write(6,*) "XYZ file(s) to be converted"
read(5,"(A)") kkchar
write(6,*) trim(kkchar)
i=0
nfich=0
nat=0
do while (.true.)
 nfich=nfich+1
 read(kkchar,*,iostat=i) (kk,j=1,nfich)
 if (i.ne.0) exit
 open(1,file=kk,status="old",iostat=j)
 if (j.ne.0) then
  stop "File does not exist"
 endif
 read(1,*) j
 close(1)
 write(6,*) "File ready ",trim(kk),j
 nat=nat+j
enddo
nfich=nfich-1

allocate(at(nat),geom(nat,3),noat(nat),mass(nat),vel(nat,3),chat(nat))

k=0
do n=1,nfich
 read(kkchar,*) (kk,j=1,n)
 open(1,file=kk,status="old")
 read(1,*) l
 read(1,*)
 do i=1,l
  read(1,"(A)") kk
  kk=trim(kk)//" 0. 0. 0. 0."
! write(6,"(A)") kkchar
  read(kk,*) at(i+k),(geom(i+k,j),j=1,3),(vel(i+k,j),j=1,3),chat(i)
  do j=1,3
   geom(i+k,j)=geom(i+k,j)/bohr2ang
   vel(i+k,j)=vel(i+k,j)/1e3
  enddo
  call getnoatandmass(at(i+k),noat(i+k),mass(i+k))
  if (noat(i+k).gt.1e-6) chat(i+k)=noat(i+k)
 enddo
 close(1)
 k=k+l
enddo

write(6,*) "Select type of output file"
read(5,*) kkchar
select case(kkchar)
 case("geom")
  write(6,*) "Creating file geom"
  open(1,file="geom")
  do i=1,nat
   write(1,"(A2,x,F5.1,4(x,F30.20))") at(i),noat(i),geom(i,:),mass(i)
  enddo
  write(1,*) ""
  do i=1,nat
   write(1,"(3(x,E30.20e3))") vel(i,:)
  enddo
 case("bagel")
  write(6,*) "Creating file geom.json"
  open(1,file="geom.json")
  write(1,"(A)") '{ "title" : "molecule", '
  write(1,"(A)") '"basis" : "",'
  write(1,"(A)") '"df_basis" : "",'
  write(1,"(A)") '"angstrom" : false,'
  write(1,"(A)") '"geometry" : [ '
  do i=1,nat
   if (i.ne.nat) then
    write(1,"(A,A,A,5(x,F20.10,A))") &
     '{"atom" : "',trim(at(i)),'","xyz" : [ ',&
      geom(i,1),', ',geom(i,2),', ',geom(i,3),&
     ' ], "mass" : ',mass(i), &
     ', "charge" : ',chat(i), &
     ' },'
   else
    write(1,"(A,A,A,4(x,F20.10,A))") &
     '{"atom" : "',trim(at(i)),'","xyz" : [ ',&
      geom(i,1),', ',geom(i,2),', ',geom(i,3),&
     ' ], "mass" : ',mass(i), &
     ', "charge" : ',chat(i), &
     ' }] }'
   endif
  enddo
 case("bagel0")
  write(6,*) "Creating file geom.json"
  open(1,file="geom.json")
  write(1,"(A)") '"geometry" : [ '
  do i=1,nat
   if (i.ne.nat) then
    write(1,"(A,A,A,5(x,F20.10,A))") &
     '{"atom" : "',trim(at(i)),'","xyz" : [ ',&
      geom(i,1),', ',geom(i,2),', ',geom(i,3),&
     ' ], "mass" : ',mass(i), &
     ' , "charge" : ',chat(i), &
     ' },'
   else
    write(1,"(A,A,A,5(x,F20.10,A))") &
     '{"atom" : "',trim(at(i)),'","xyz" : [ ',&
      geom(i,1),', ',geom(i,2),', ',geom(i,3),&
     ' ], "mass" : ',mass(i), &
     ' , "charge" : ',chat(i), &
     ' }]'
   endif
  enddo
 case("molden")
  write(6,*) "Creating file geom.molden"
  open(1,file="geom.molden")
  write(1,"(A)") '[Molden format] '
  write(1,"(A)") '[Atoms] (AU)'
  do i=1,nat
   write(1,"(A,2(X,I0),3(x,F20.10))") &
    trim(at(i)),i,int(noat(i)),&
     geom(i,1),geom(i,2),geom(i,3)
  enddo
 case("gamess")
  write(6,*) "Creating file geom.gamess (nosym)"
  open(1,file="geom.gamess")
  write(1,"(A)") " $DATA"
  write(1,*) ""
  write(1,*) "C1"
  do i=1,nat
   write(1,"(A,X,I0,3(x,F20.10))") & 
   trim(at(i)),int(noat(i)),&
    (geom(i,j)*bohr2ang,j=1,3)
  enddo
  write(1,*) "$END"
  close(1)
 case default
  write(6,*) "Valid output files are geom, bagel, bagel0, molden or gamess"
end select

end

subroutine getnoatandmass(at,noat,mass)
implicit none
character*2, intent(inout) :: at
real*8, intent(out) :: noat,mass

integer i
character*2 sat(0:118),sat2(0:118)
real*8 atmass(0:118)

sat(0)="XX"
sat(1)="H"
sat(2)="He"
sat(3)="Li"
sat(4)="Be"
sat(5)="B"
sat(6)="C"
sat(7)="N"
sat(8)="O"
sat(9)="F"
sat(10)="Ne"
sat(11)="Na"
sat(12)="Mg"
sat(13)="Al"
sat(14)="Si"
sat(15)="P"
sat(16)="S"
sat(17)="Cl"
sat(18)="Ar"
sat(19)="K"
sat(20)="Ca"
sat(21)="Sc"
sat(22)="Ti"
sat(23)="V"
sat(24)="Cr"
sat(25)="Mn"
sat(26)="Fe"
sat(27)="Co"
sat(28)="Ni"
sat(29)="Cu"
sat(30)="Zn"
sat(31)="Ga"
sat(32)="Ge"
sat(33)="As"
sat(34)="Se"
sat(35)="Br"
sat(36)="Kr"
sat(37)="Rb"
sat(38)="Sr"
sat(39)="Y"
sat(40)="Zr"
sat(41)="Nb"
sat(42)="Mo"
sat(43)="Tc"
sat(44)="Ru"
sat(45)="Rh"
sat(46)="Pd"
sat(47)="Ag"
sat(48)="Cd"
sat(49)="In"
sat(50)="Sn"
sat(51)="Sb"
sat(52)="Te"
sat(53)="I"
sat(54)="Xe"
sat(55)="Cs"
sat(56)="Ba"
sat(57)="La"
sat(58)="Ce"
sat(59)="Pr"
sat(60)="Nd"
sat(61)="Pm"
sat(62)="Sm"
sat(63)="Eu"
sat(64)="Gd"
sat(65)="Tb"
sat(66)="Dy"
sat(67)="Ho"
sat(68)="Er"
sat(69)="Tm"
sat(70)="Yb"
sat(71)="Lu"
sat(72)="Hf"
sat(73)="Ta"
sat(74)="W"
sat(75)="Re"
sat(76)="Os"
sat(77)="Ir"
sat(78)="Pt"
sat(79)="Au"
sat(80)="Hg"
sat(81)="Tl"
sat(82)="Pb"
sat(83)="Bi"
sat(84)="Po"
sat(85)="At"
sat(86)="Rn"
sat(87)="Fr"
sat(88)="Ra"
sat(89)="Ac"
sat(90)="Th"
sat(91)="Pa"
sat(92)="U"
sat(93)="Np"
sat(94)="Pu"
sat(95)="Am"
sat(96)="Cm"
sat(97)="Bk"
sat(98)="Cf"
sat(99)="Es"
sat(100)="Fm"
sat(101)="Md"
sat(102)="No"
sat(103)="Lr"
sat(104)="Rf"
sat(105)="Db"
sat(106)="Sg"
sat(107)="Bh"
sat(108)="Hs"
sat(109)="Mt"
sat(110)="Ds"
sat(111)="Rg"
sat(112)="Cn"
sat(113)="Nh"
sat(114)="Fl"
sat(115)="Mc"
sat(116)="Lv"
sat(117)="Ts"
sat(118)="Og"
sat2(0)="xx"
sat2(1)="h"
sat2(2)="he"
sat2(3)="li"
sat2(4)="be"
sat2(5)="b"
sat2(6)="c"
sat2(7)="n"
sat2(8)="o"
sat2(9)="f"
sat2(10)="ne"
sat2(11)="na"
sat2(12)="mg"
sat2(13)="al"
sat2(14)="si"
sat2(15)="p"
sat2(16)="s"
sat2(17)="cl"
sat2(18)="ar"
sat2(19)="k"
sat2(20)="ca"
sat2(21)="sc"
sat2(22)="ti"
sat2(23)="v"
sat2(24)="cr"
sat2(25)="mn"
sat2(26)="fe"
sat2(27)="co"
sat2(28)="ni"
sat2(29)="cu"
sat2(30)="zn"
sat2(31)="ga"
sat2(32)="ge"
sat2(33)="as"
sat2(34)="se"
sat2(35)="br"
sat2(36)="kr"
sat2(37)="rb"
sat2(38)="sr"
sat2(39)="y"
sat2(40)="zr"
sat2(41)="nb"
sat2(42)="mo"
sat2(43)="tc"
sat2(44)="ru"
sat2(45)="rh"
sat2(46)="pd"
sat2(47)="ag"
sat2(48)="cd"
sat2(49)="in"
sat2(50)="sn"
sat2(51)="sb"
sat2(52)="te"
sat2(53)="i"
sat2(54)="xe"
sat2(55)="cs"
sat2(56)="ba"
sat2(57)="la"
sat2(58)="ce"
sat2(59)="pr"
sat2(60)="nd"
sat2(61)="pm"
sat2(62)="sm"
sat2(63)="eu"
sat2(64)="gd"
sat2(65)="tb"
sat2(66)="dy"
sat2(67)="ho"
sat2(68)="er"
sat2(69)="tm"
sat2(70)="yb"
sat2(71)="lu"
sat2(72)="hf"
sat2(73)="ta"
sat2(74)="w"
sat2(75)="re"
sat2(76)="os"
sat2(77)="ir"
sat2(78)="pt"
sat2(79)="au"
sat2(80)="hg"
sat2(81)="tl"
sat2(82)="pb"
sat2(83)="ni"
sat2(84)="po"
sat2(85)="at"
sat2(86)="rn"
sat2(87)="fr"
sat2(88)="ra"
sat2(89)="ac"
sat2(90)="th"
sat2(91)="pa"
sat2(92)="u"
sat2(93)="np"
sat2(94)="pu"
sat2(95)="am"
sat2(96)="cm"
sat2(97)="bk"
sat2(98)="cf"
sat2(99)="es"
sat2(100)="fm"
sat2(101)="md"
sat2(102)="no"
sat2(103)="lr"
sat2(104)="rf"
sat2(105)="db"
sat2(106)="sg"
sat2(107)="bh"
sat2(108)="hs"
sat2(109)="mt"
sat2(110)="ds"
sat2(111)="rg"
sat2(112)="cn"
sat2(113)="nh"
sat2(114)="fl"
sat2(115)="mc"
sat2(116)="lv"
sat2(117)="ts"
sat2(118)="og"
atmass(0)=0.
atmass(1)=1.00782503223
atmass(2)=3.0160293201
atmass(3)=6.0151228874
atmass(4)=9.012183065
atmass(5)=10.01293695
atmass(6)=12.0000000
atmass(7)=14.00307400443
atmass(8)=15.99491461957
atmass(9)=18.99840316273
atmass(10)=19.9924401762
atmass(11)=22.9897692820
atmass(12)=23.985041697
atmass(13)=26.98153853
atmass(14)=27.97692653465
atmass(15)=30.97376199842
atmass(16)=31.9720711744
atmass(17)=34.968852682
atmass(18)=35.967545105
atmass(19)=38.9637064864
atmass(20)=39.962590863
atmass(21)=44.95590828
atmass(22)=45.95262772
atmass(23)=49.94715601
atmass(24)=49.94604183
atmass(25)=54.93804391
atmass(26)=53.93960899
atmass(27)=58.93319429
atmass(28)=57.93534241
atmass(29)=62.92959772
atmass(30)=63.92914201
atmass(31)=68.9255735
atmass(32)=69.92424875
atmass(33)=74.92159457
atmass(34)=73.922475934
atmass(35)=78.9183376
atmass(36)=77.92036494
atmass(37)=84.9117897379
atmass(38)=83.9134191
atmass(39)=88.9058403
atmass(40)=89.9046977
atmass(41)=92.9063730
atmass(42)=91.90680796
atmass(43)=96.9063667
atmass(44)=95.90759025
atmass(45)=102.9054980
atmass(46)=101.9056022
atmass(47)=106.9050916
atmass(48)=105.9064599
atmass(49)=112.90406184
atmass(50)=111.90482387
atmass(51)=120.9038120
atmass(52)=119.9040593
atmass(53)=126.9044719
atmass(54)=123.9058920
atmass(55)=132.9054519610
atmass(56)=129.9063207
atmass(57)=137.9071149
atmass(58)=135.90712921
atmass(59)=140.9076576
atmass(60)=141.9077290
atmass(61)=144.9127559
atmass(62)=143.9120065
atmass(63)=150.9198578
atmass(64)=151.9197995
atmass(65)=158.9253547
atmass(66)=155.9242847
atmass(67)=164.9303288
atmass(68)=161.9287884
atmass(69)=168.9342179
atmass(70)=167.9338896
atmass(71)=174.9407752
atmass(72)=173.9400461
atmass(73)=179.9474648
atmass(74)=179.9467108
atmass(75)=184.9529545
atmass(76)=183.9524885
atmass(77)=190.9605893
atmass(78)=189.9599297
atmass(79)=196.96656879
atmass(80)=195.9658326
atmass(81)=202.9723446
atmass(82)=203.9730440
atmass(83)=208.9803991
atmass(84)=208.9824308
atmass(85)=209.9871479
atmass(86)=210.9906011
atmass(87)=223.0197360
atmass(88)=223.0185023
atmass(89)=227.0277523
atmass(90)=230.0331341
atmass(91)=231.0358842
atmass(92)=233.0396355
atmass(93)=236.046570
atmass(94)=238.0495601
atmass(95)=241.0568293
atmass(96)=243.0613893
atmass(97)=247.0703073
atmass(98)=249.0748539
atmass(99)=252.082980
atmass(100)=257.0951061
atmass(101)=258.0984315
atmass(102)=259.10103
atmass(103)=262.10961
atmass(104)=267.12179
atmass(105)=268.12567
atmass(106)=271.13393
atmass(107)=272.13826
atmass(108)=270.13429
atmass(109)=276.15159
atmass(110)=281.16451
atmass(111)=280.16514
atmass(112)=285.17712
atmass(113)=284.17873
atmass(114)=289.19042
atmass(115)=288.19274
atmass(116)=293.20449
atmass(117)=292.20746
atmass(118)=294.21392

do i=1,118
 if (at==sat(i) .or. at==sat2(i)) then
  noat=float(i)
  mass=atmass(i)
  at=sat(i)
 endif
enddo


end

