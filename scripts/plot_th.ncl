load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;=================================
begin
;=================================
; Settings.
;---------------------------------

;nx = 301
;nz = 80
;nt = 73

;nx = 801
nx = 201
;nz = 80
nz = 31
;nt = 73
;nt = 37

;nx = 501
;nz = 200
nt = 180

;nt = 1
;nt = 2
;nt = 8
;nt = 15
;nt = 28
;nt = 35
;nt = 40
;nt = 61
;nt = 73
;nt = 91
;nt = 121
;nt = 181
;nt = 141
;nt = 301
;nt = 401
;nt = 1001
;var = "theta"
;var = "theta_M_1"
;var = "cldfra"
;var = "theta_1"
;var = "theta_M"
;var = "pi_1"
;var = "u"
var = "w"
;var = "qv"
;var = "qc"
;var = "qr"
;var = "qi"
;var = "qs"
;var = "qg"
halo = 4
;=================================
; Read data.
;---------------------------------
setfileoption ("bin", "ReadByteOrder", "BigEndian")
;data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz,nx/),"double")
data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz+2*halo,nx+2*halo/),"double")
;data = data*1000.
;data = smth9(data,.5,.25,True)
;print(data(nt-1,10,:))
;=================================
; Plot.
;---------------------------------
wks = gsn_open_wks ("pdf","th")                  ; open workstation
;wks = gsn_open_wks ("pdf","tb_w")                  ; open workstation
;wks = gsn_open_wks ("pdf","tb_u")                  ; open workstation
gsn_define_colormap(wks,"temp_diff_18lev")    ; For theta
;gsn_define_colormap(wks,"BlueWhiteOrangeRed")  ; For u & w

res = True

res@gsnMaximize = True
res@vpWidthF  = 1.0                    ; change the aspect ratio
;res@vpHeightF = 0.1
res@vpHeightF = 0.28
;res@vpWidthF  = 0.5                    ; change the aspect ratio
;res@vpHeightF = 1.0
;res@vpXF      = .1                     ; location of where plot starts
;res@vpYF      = .7
res@tmXBLabelFontHeightF = 0.015       ; label font height
res@tmYLLabelFontHeightF = 0.015

res@tmXTOn = False
res@tmYROn = False

;res@cnLevelSelectionMode           = "ExplicitLevels"
;res@cnLevels           = (/0.01,0.02,0.04,0.08,0.16,0.32,0.64,1.28,2.56,5.12,10.24/)
res@cnLineLabelsOn = False
res@cnInfoLabelOn = True
;res@cnInfoLabelOn = False

res@tiMainString = "Thunderstorm"
;res@tiXAxisString      = "X"  
;res@tiYAxisString      = "Z (height)"
res@gsnContourNegLineDashPattern = 2

res@tiMainFontHeightF  = 0.015                    ; change font heights
;res@tiXAxisFontHeightF  = 0.015                    ; change font heights
;res@tiYAxisFontHeightF  = 0.015                    ; change font heights

;res@gsnLeftString   = "dx = 50m; dz = 50m; dt = 0.1s"               ; add the gsn titles
res@gsnLeftString   = "dx = 250m; dz = 250m; dt = 0.3s"               ; add the gsn titles
res@gsnLeftStringFontHeightF   = 0.015

;plot = gsn_contour(wks,pi_data,res)
plots = new(nt,graphic)
;do i = 0, nt-1, 100
do i = 0, nt-1
	print("Step: "+i+" plotting...")
	;res@gsnCenterString = i
	;res@gsnCenterStringFontHeightF   = 0.015
	k = i*100
	;k = i*1000
	res@gsnRightString  = "Step: " + k
	;res@gsnRightString  = "Step: " + i
	res@gsnRightStringFontHeightF   = 0.015
	;plots(0) = gsn_csm_contour(wks,data(i,:,:),res)
	plots(0) = gsn_csm_contour(wks,data(i,halo:halo+nz-1,halo:halo+nx-1),res)
end do

;=================================
end
;=================================
