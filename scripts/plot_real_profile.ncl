load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;=================================
begin
;=================================
; Settings.
;---------------------------------

nx = 31
nz = 31
nt = 10

;var = "theta"
;var = "theta_M_1"
;var = "cldfra"
;var = "theta_1"
;var = "theta_M"
;var = "pi_1"
;var = "w"
;var = "qv"
;var = "qc"
;var = "qr"
;var = "qi"
;var = "qs"
;var = "qg"
halo = 4
;rain = cbinread("../output/modelvar_rain.bin",(/nt,nx+2*halo/),"double")
;printVarSummary(rain)
;print(rain(nt-1,:))
nlat = 100
nlon = 100
nlev = 32
;=================================
; Read data.
;---------------------------------
setfileoption ("bin", "ReadByteOrder", "BigEndian")

zh_full = cbinread("../output/zh_full.dat",(/2,nlev,nlat,nlon/),"float")
zh_half = cbinread("../output/zh_half.dat",(/2,nlev-1,nlat,nlon/),"float")

tpi = 50 - 1
tpj = 50 - 1
time = 1 - 1
zz_full = zh_full(time,:,tpj,tpi)
zz_half = zh_half(time,:,tpj,tpi)

time = nt - 1

var = "u"
data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz+2*halo,nx+2*halo/),"double")
u_ave = dim_avg(data)
uu = u_ave(time,halo:nz+halo-1)
delete(data)
;printVarSummary(u_ave)

var = "w"
data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz+2*halo,nx+2*halo/),"double")
w_ave = dim_avg(data)
ww = w_ave(time,halo:nz+halo)
delete(data)

var = "theta"
data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz+2*halo,nx+2*halo/),"double")
theta_ave = dim_avg(data)
tt = theta_ave(time,halo:nz+halo)
delete(data)

var = "qv"
data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz+2*halo,nx+2*halo/),"double")
qv_ave = dim_avg(data)
qqv = qv_ave(time,halo:nz+halo)
delete(data)

var = "qc"
data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz+2*halo,nx+2*halo/),"double")
qc_ave = dim_avg(data)
qqc = qc_ave(time,halo:nz+halo)
delete(data)

var = "qr"
data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz+2*halo,nx+2*halo/),"double")
qr_ave = dim_avg(data)
qqr = qr_ave(time,halo:nz+halo)
delete(data)

var = "qi"
data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz+2*halo,nx+2*halo/),"double")
qi_ave = dim_avg(data)
qqi = qi_ave(time,halo:nz+halo)
delete(data)

var = "qs"
data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz+2*halo,nx+2*halo/),"double")
qs_ave = dim_avg(data)
qqs = qs_ave(time,halo:nz+halo)
delete(data)

var = "qg"
data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz+2*halo,nx+2*halo/),"double")
qg_ave = dim_avg(data)
qqg = qg_ave(time,halo:nz+halo)
delete(data)
;=================================
; Plot.
;---------------------------------
wks = gsn_open_wks ("pdf","real")                  ; open workstation

profile = True
profile@gsnMaximize = True
profile@vpWidthF  = 0.4                    ; change the aspect ratio
profile@vpHeightF = 1.0

;profile@tiXAxisFontHeightF  = 0.025                    ; change font heights
;profile@tiYAxisFontHeightF  = 0.025                    ; change font heights

;;profile@tmXBMode        = "Explicit"
;;profile@tmXBValues      = ispan(0,300,50)
;;profile@tmXBLabels      = (/"0","50","100","150","200","250","300"/)
;;profile@tmXBLabelFontHeightF = 0.02
;;profile@tmXBMinorValues = ispan(0,300,10)

;profile@xyLineThicknessF = 2.
;profile@tmYLLabelFontHeightF   = 0.020
;profile@tmXBLabelFontHeightF   = 0.020

profile@trYMinF = 0
profile@trYMaxF = 30000
;;profile@gsnYRefLine = 0.
;;profile@gsnYRefLineThicknessF = 4.

profile@tmXTOn = False
profile@tmYROn = False

profile@gsnDraw = False
profile@gsnFrame = False

plots = new(6,graphic)

;profile@tiXAxisString      = "Potential temperature (K)"
profile@tiYAxisString      = "Height (km)"  
i = 0
;plots(i) = gsn_csm_xy(wks,theta,zz_full,profile)
;i = i+1
;print("OK theta")
;delete(profile@tiYAxisString)

profile@tiXAxisString      = "u (m s-1)"
profile@gsnXRefLine = 0.
profile@gsnXRefLineThicknessF = 2.
plots(i) = gsn_csm_xy(wks,uu,zz_half,profile)
i = i+1
print("OK u")

delete(profile@tiYAxisString)
profile@tiXAxisString      = "w (m s-1)"
profile@gsnXRefLine = 0.
profile@gsnXRefLineThicknessF = 2.
plots(i) = gsn_csm_xy(wks,ww,zz_full,profile)
i = i+1
print("OK w")

;profile@tiXAxisString      = "Water vapor (g kg-1)"
profile@tiXAxisString      = "Potential temperature (K)"
plots(i) = gsn_csm_xy(wks,tt,zz_full,profile)
i = i+1
print("OK pi")


profile@tiYAxisString      = "Height (km)"  
profile@tiXAxisString      = "qv (g kg-1)"
plots(i) = gsn_csm_xy(wks,qqv,zz_full,profile)
i = i+1
print("OK qv")

;=============================================
profile@tiXAxisString      = "qc, qi, qr, qs, and qg (g kg-1)"
profile@xyDashPatterns     = (/0,2,0,4,1/)
profile@xyLineThicknesses = (/2.,1.,1.,1.,1.,1.,4./)
profile@xyExplicitLegendLabels = (/" qg"," qs"," qr"," qi", " qc"/)

profile@lgPerimOn              = True
profile@pmLegendDisplayMode    = "Always"   ; Display a legend.
profile@pmLegendWidthF         = 0.3        ; Make it smaller
profile@pmLegendHeightF        = 0.3        ; Make it smaller
profile@pmLegendOrthogonalPosF = -1.00      ; Move into to plot
profile@pmLegendParallelPosF   = 0.55      ; Move to right
profile@lgLabelFontHeightF     = 0.020

qq = (/qqg,qqs,qqr,qqi,qqc/)
plots(i) = gsn_csm_xy(wks,qq,zz_full,profile)
print("OK qq")
;=============================================

;;delete(profile@tiYAxisString)
;;profile@tiXAxisString      = "qc (g kg-1)"
;;plots(i) = gsn_csm_xy(wks,qqc,zz_full,profile)
;;i = i+1
;;print("OK qc")

;;delete(profile@tiYAxisString)
;;profile@tiXAxisString      = "qr (g kg-1)"
;;plots(i) = gsn_csm_xy(wks,qqr,zz_full,profile)
;;i = i+1
;;print("OK qr")

;;delete(profile@tiYAxisString)
;;profile@tiXAxisString      = "qi (g kg-1)"
;;plots(i) = gsn_csm_xy(wks,qqi,zz_full,profile)
;;i = i+1
;;print("OK qi")

;;delete(profile@tiYAxisString)
;;profile@tiXAxisString      = "qs (g kg-1)"
;;plots(i) = gsn_csm_xy(wks,qqs,zz_full,profile)
;;i = i+1
;;print("OK qs")

;;delete(profile@tiYAxisString)
;;profile@tiXAxisString      = "qg (g kg-1)"
;;plots(i) = gsn_csm_xy(wks,qqg,zz_full,profile)
;;i = i+1
;;print("OK qg")

pres =True
pres@gsnMaximize = True
pres@gsnPanelBottom   = 0.05			 ; shrink panel plot by setting bottom edge of plot
pres@gsnPanelTop      = 1.00			 ; shrink panel plot by setting top edge of plot
pres@amJust = "TopLeft"
pres@gsnPanelFigureStrings = (/"a)","b)","c)","d)","e)","f)","g)","h)","i)","j)","k)"/)      ; add strings to panel
pres@gsnPanelFigureStringsFontHeightF = 0.015

gsn_panel(wks,plots,(/2,3/),pres)

;=================================
end
;=================================
