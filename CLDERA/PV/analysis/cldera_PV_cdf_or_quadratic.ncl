;*************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
;************************************************
begin

;=====
; declare the file name of the plot, data file
;=====

; these need to be defined as cmd args:
; res_name     -- e.g. 'ne16np4'
; compset_name -- e.g. 'FIDEAL'
; grid_name    -- 'dyn' or 'phys'
; remap_res    -- e.g. '91x180'

plotname = "pv_verify_" + compset_name + "_" +res_name
case_name= "E3SM_" + res_name + "_L72_" + compset_name + "_PV"
dir = "/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/pv_cases/"
data = dir + case_name + "/run/" + case_name + \
       ".eam.h0.0001-01-01-00000.regrid." +grid_name + "." + remap_res + "_bilinear.nc"

;=====
; declare a label for the dycore that appears in the plot
;=====
 dycore   = "HOMME"               ; label in plot and filename

;=====
; declare a label for the type of dissipation
;=====
 diff   = "PV_homme_internal_quadratic_NCL_recon_NCL_built-in"   ; 4th-order divergence damping 

;+++++++++++++++++++++++++++++++++++++++++++++++++++++
; day, lev
;+++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  pday = 30                     ; default day 30
  ;plev =  5                    ; top
  ;plev =  96.47                ; ~tropopause
  plev =  180                   ; surface
  ;plev =  950                  ; surface
  ;type = "X11"                 ; default X11 window
  type = "eps"                  ; default X11 window

;=====
; complete the plotname
;=====
  plotname = plotname+"_"+grid_name+"Grid_"+diff+"_near_"+plev+"_day"+pday+"_cfd_homme"


;************************************************
; open netCDF file
;************************************************
  a = addfile(data,"r")

;************************************************
; read variables
;************************************************
  lev  = a->lev
  lat  = a->lat
  klev  = dimsizes(lev) - 1              ; index counts from 0 ... klev
  hyam = a->hyam
  hybm = a->hybm
  hyai = a->hyai
  hybi = a->hybi

  p0   = 100000.d0      ; Pa
  p0mb = p0*0.01d0      ; hPa
  gravity = 9.80616d0   ; gravitational acceleration in m/s^2
  sidereal_day = 86164.d0
  exponent = 2.d0/7.d0
  pi = 4.d0*atan(1.d0)
  deg_rad = pi/180.d0
  omega = 2.d0*pi/sidereal_day
  rearth = 6371000.d0   ; radius in m

;*************
; read data
;*************
  ps = a->PS(pday,:,:)
  ps@units      = "Pa"
  ps@long_name  = "Surface pressure"

  t = a->T(pday,:,:,:)
  t@units = "K"
  u = a->U(pday,:,:,:)
  v = a->V(pday,:,:,:)

  pv = a->PV(pday,:,:,:)
  pv = pv*1.d6                ; regular PV scaling
  pv@long_name = "Potential vorticity (internal HOMME, quadratic)"
  pv@units     = "PVU"
  printMinMax(pv({plev},:,:),True)
  pv_min = tofloat(min(pv({plev},:,:)))
  pv_max = tofloat(max(pv({plev},:,:)))

; create constant 2D pressure field
  p00 = ps
  p00(:,:) = p0

; compute eta
  eta = hyam
  eta = hyam + hybm
  eta@long_name  = "Eta levels"
  printMinMax(eta,True)

; pressure at full model level
  pm = pres_hybrid_ccm(ps,p0,hyam,hybm) 
  copy_VarCoords(u,pm)
  pm@units      = "Pa"
  pm@long_name  = "Pressure"
  printVarSummary (pm)
  printMinMax(pm({plev},:,:),True)
  
; potential temperature
  theta = t
  theta = t*(p0/pm)^exponent 
  theta@long_name  = "Potential temperature"
  theta@units  = "K"
  printMinMax(theta({plev},:,:),True)
  
; horizontal gradients of the potential temperature
  grad_lon=t
  grad_lat=t
  gradsf (theta, grad_lon, grad_lat)

; finite difference computation is an alternative
;  gradLatLon  = grad_latlon_cfd (theta, theta&lat, theta&lon, True, False)
;
;  grad_lat  = gradLatLon[0]    ; for clarity; explicitly extract variables from returned 'list'
;  grad_lon  = gradLatLon[1]
;  delete(gradLatLon)

  grad_lon@long_name = "longitudinal gradient (derivative)"
  grad_lat@long_name = "latitudinal gradient (derivative)"
  grad_lat@units     = "K/m"
  grad_lon@units     = "K/m"

; Potential Vorticity via NCL routine
  gridType = 1  ; Equidistant grid
  pv_NCL_internal = pot_vort_hybrid(pm,u,v,t,lat, gridType,  0)
  pv_NCL_internal = pv_NCL_internal * 1.d6
  pv_NCL_internal@long_name = "Potential vorticity (built-in NCL)"
  pv_NCL_internal@units     = "PVU"
  printMinMax(pv_NCL_internal({plev},:,:),True)
  pv_NCL_internal_min = tofloat(min(pv_NCL_internal({plev},:,:)))
  pv_NCL_internal_max = tofloat(max(pv_NCL_internal({plev},:,:)))


;------------
;  compute relative vorticity 
;------------
  vort = uv2vrF_Wrap (u,v)
  copy_VarMeta(u,vort)
  vort@long_name = "Relative vorticity (built-in NCL)"
  vort@units     = "10~S~-5~N~ 1/s"

; compute the Coriolis parameter
  f = lat
  f = 2.d0*omega*sin(lat*deg_rad)
  ff = conform_dims(dimsizes(u), f, 1)
 
; compute the PV via centered finite differences in the vertical, omit the top and bottom level
  pv_NCL_external = pv_NCL_internal
  pv_NCL_external@long_name = "Potential vorticity (NCL)"

  do k = 1, klev-1

; correct equation
    pv_NCL_external(k,:,:) = - gravity/((hyam(k+1)-hyam(k-1))*p00 + (hybm(k+1)-hybm(k-1))*ps) * ((vort(k,:,:)+ff(k,:,:))*(theta(k+1,:,:)-theta(k-1,:,:)) - (v(k+1,:,:)-v(k-1,:,:))*grad_lon(k,:,:) + (u(k+1,:,:)-u(k-1,:,:))*grad_lat(k,:,:))  

  end do

; scale
  pv_NCL_external = pv_NCL_external*1.d6

; compute difference: HOMME - NCL PV
  pv_difference = pv
  pv_difference = pv - pv_NCL_external
  pv_difference@long_name = "PV difference (HOMME - NCL)"
  pv_difference_min = tofloat(min(pv_difference({plev},:,:)))
  pv_difference_max = tofloat(max(pv_difference({plev},:,:)))

  vort            = vort*1.d5
  printMinMax(pv_NCL_external({plev},:,:),True)
  pv_NCL_external_min = tofloat(min(pv_NCL_external({plev},:,:)))
  pv_NCL_external_max = tofloat(max(pv_NCL_external({plev},:,:)))

;************************************************
;plot resources [options]
;************************************************
  plot = new(4,graphic)

  pres = True               ; panel
  pres@txString    = case_name + ", " + grid_name + " grid, day "+pday+", model level near "+plev+" hPa"
  pres@gsnPanelYWhiteSpacePercent = 0
  pres@gsnPanelXWhiteSpacePercent = 3
  pres@gsnFrame         = False
  pres@gsnPanelFigureStringsFontHeightF = 0.012
  pres@amJust           = "TopLeft"
  pres@gsnPanelTop = 0.95
  pres@txFontHeightF = 0.014

  res = True  
  res@gsnDraw           = False              ; panel ... set to False
  res@gsnFrame          = False
  res@cnFillOn          = True
  res@gsnSpreadColors   = True
  res@lbLabelAutoStride = False
  res@tiMainString      = ""
  res@vpWidthF          = 0.44
  res@vpHeightF         = 0.22
  res@cnLinesOn         = False
  res@gsnCenterString   = " "
  res@cnLevelSelectionMode = "ManualLevels"
  res@pmLabelBarWidthF = 0.43
  res@gsnSpreadColors   = True
  res@gsnStringFontHeightF = 0.012
  res@gsnSpreadColorStart = 5

;=======
; open plot and define color table
;=======
  wks = gsn_open_wks(type,plotname)   
; gsn_define_colormap(wks,"WhBlGrYeRe")          ; select color table
  gsn_define_colormap(wks,"amwg_blueyellowred")  ; select color table
; gsn_define_colormap(wks,"cmp_flux")            ; select color table
; gsn_define_colormap(wks,"prcp_2")              ; select color table


  res@cnLevelSelectionMode = "AutomaticLevels"
  res@lbLabelAngleF    =    90
  ;res@cnLevelSelectionMode = "ManualLevels"
  ;res@cnLevelSpacingF = 0.3
  ;res@cnMinLevelValF  = -1.2
  ;res@cnMaxLevelValF  = 1.2
  plot(2) = gsn_csm_contour(wks,pv_difference({plev},:,:),res)


  res@cnLevelSelectionMode = "AutomaticLevels"
  res@lbLabelAngleF    =    90
  ;res@cnLevelSelectionMode = "ExplicitLevels"
  ;res@cnLevels = (/ -1000., -300., -100., -30., -10., -3., -1., 0., \
  ;                   1., 3., 10., 30., 100., 300., 1000./)
  ;res@cnLevels = (/ -20., -15., -10., -7., -5., -1., 0., \
  ;                   1., 5., 7., 10., 15., 20. /)
  plot(0) = gsn_csm_contour(wks,pv({plev},:,:),res)
  plot(1) = gsn_csm_contour(wks,pv_NCL_internal({plev},:,:),res)
  plot(3) = gsn_csm_contour(wks,pv_NCL_external({plev},:,:),res)

  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.01                     ; font smaller. default big
  gsn_text_ndc(wks,"min: " + pv_min + ", max: " + pv_max,0.25,0.47,txres) 
  gsn_text_ndc(wks,"min: "+pv_NCL_internal_min + ", max: " + pv_NCL_internal_max,0.75,0.47,txres) 
  gsn_text_ndc(wks,"min: "+pv_NCL_external_min + ", max: " + pv_NCL_external_max,0.75,0.11,txres) 
  gsn_text_ndc(wks,"min: "+pv_difference_min + ", max: " + pv_difference_max,0.25,0.11,txres) 

  gsn_panel(wks,plot,(/2,2/),pres)
  frame(wks)

  ;trim the white space
  ;if (type .ne. "X11") then
  ;  system("convert -trim "+plotname+"."+type+" "+plotname+"."+type)
  ;end if

end
