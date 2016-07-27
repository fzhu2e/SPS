load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;=================================
begin
;=================================
; Settings.
;---------------------------------

;nx = 100
;nz = 200

;nx = 10
;nz = 20

nx = 401
nz = 200
;nt = 1
;nt = 3
;nt = 5
;nt = 1
;nt = 11
;nt = 25
;nt = 55
;nt = 75
;nt = 91
nt = 101
;nt = 181
;nt = 141
;nt = 301
;nt = 401
;nt = 1001
;var = "theta"
;var = "theta_M"
var = "theta_M_1"
;var = "pi_1"
;var = "u"
;var = "w"
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
;data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz,nx/),"float")
data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz+2*halo,nx+2*halo/),"double")
;print(data)
;=================================
; Plot.
;---------------------------------
wks = gsn_open_wks ("pdf","wb")                  ; open workstation
gsn_define_colormap(wks,"temp_diff_18lev")    ; For theta
;gsn_define_colormap(wks,"BlueWhiteOrangeRed")  ; For u & w

res = True

res@gsnMaximize = True
res@vpWidthF  = 0.6                    ; change the aspect ratio
res@vpHeightF = 1.0
;res@vpWidthF  = 0.5                    ; change the aspect ratio
;res@vpHeightF = 1.0
;res@vpXF      = .1                     ; location of where plot starts
;res@vpYF      = .7
res@tmXBLabelFontHeightF = 0.015       ; label font height
res@tmYLLabelFontHeightF = 0.015

res@tmXTOn = False
res@tmYROn = False


; LabelBar
res@cnFillOn = True
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
;res@pmLabelBarDisplayMode   = "Always" ; turn on color bar
res@pmLabelBarSide          = "Bottom" ; put it on the bottom of plot
res@pmLabelBarWidthF        = 1.0      ; control size of colorbar
res@pmLabelBarHeightF       = 0.1      ;
;res@pmLabelBarOrthogonalPosF = -0.02   ; position wrt plot
res@lbLabelFontHeightF      = 0.015     ; label bar font
res@lbOrientation           = "Horizontal" ; label bar orientation
res@lbPerimOn               = False    ; no box around label bar
res@lbAutoManage            = False	 ; we control how drawn not plot 
;res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
;res@cnLevels            = (//)

res@gsnSpreadColors = True
res@gsnSpreadColorStart = 11    ; For theta
res@gsnSpreadColorEnd = 60

;res@cnLevelSelectionMode = "ManualLevels"   ; set explicit contour levels
; For theta
;res@cnMinLevelValF          = 0.1
;res@cnMaxLevelValF          = 1.0
;res@cnLevelSpacingF         = 0.1

; For u & w
;res@cnMinLevelValF          = - 10
;res@cnMaxLevelValF          = 10
;res@cnLevelSpacingF         = 1

;res@cnFillColors = ispan(250,10,16)

; Line
;res@cnFillOn = False
;res@cnLinesOn = True
;res@cnLineLabelsOn = True
;res@cnMonoLevelFlag = False
;res@cnLevelFlags = "LineAndLabel"

res@tiMainString = "Wet Bubble"
;res@tiXAxisString      = "X"  
;res@tiYAxisString      = "Z (height)"

res@tiMainFontHeightF  = 0.015                    ; change font heights
;res@tiXAxisFontHeightF  = 0.015                    ; change font heights
;res@tiYAxisFontHeightF  = 0.015                    ; change font heights

;res@gsnLeftString   = "dx = 50m; dz = 50m; dt = 0.1s"               ; add the gsn titles
res@gsnLeftString   = "dx = 5m; dz = 5m; dt = 0.01s"               ; add the gsn titles
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
