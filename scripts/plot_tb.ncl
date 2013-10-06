load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;=================================
begin
;=================================
; Settings.
;---------------------------------

nx = 200
nz = 100

;nx = 400
;nz = 200
;nt = 1
;nt = 3
;nt = 2
;nt = 1
;nt = 11
;nt = 65
;nt = 91
nt = 101
;nt = 141
;nt = 401
;nt = 1001
var = "theta_1"
;var = "pi_1"
;var = "u"
;var = "w"
;=================================
; Read data.
;---------------------------------
setfileoption ("bin", "ReadByteOrder", "BigEndian")
;data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz,nx/),"float")
data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz,nx/),"double")
;print(data)
pi_data = (data(:,0:nz-2,:) + data(:,1:nz-1,:))/2.
;=================================
; Plot.
;---------------------------------
wks = gsn_open_wks ("pdf","tb")                  ; open workstation
gsn_define_colormap(wks,"temp_diff_18lev")
;gsn_define_colormap(wks,"WhiteYellowOrangeRed")
;gsn_define_colormap(wks,"WhiteBlue")

res = True

res@gsnMaximize = True
res@vpWidthF  = 1.0                    ; change the aspect ratio
res@vpHeightF = 0.5
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
res@gsnSpreadColorStart = 11
res@gsnSpreadColorEnd = 60

res@cnLevelSelectionMode = "ManualLevels"   ; set explicit contour levels
res@cnMinLevelValF          = 0.2
res@cnMaxLevelValF          = 1.8
res@cnLevelSpacingF         = 0.2
;res@cnMinLevelValF          = 0.10
;res@cnMaxLevelValF          = 0.20
;res@cnLevelSpacingF         = 0.001

;res@cnFillColors = ispan(250,10,16)

; Line
;res@cnFillOn = False
;res@cnLinesOn = True
;res@cnLineLabelsOn = True
;res@cnMonoLevelFlag = False
;res@cnLevelFlags = "LineAndLabel"

res@tiMainString = "Thermal Bubble"
;res@tiXAxisString      = "X"  
;res@tiYAxisString      = "Z (height)"

res@tiMainFontHeightF  = 0.015                    ; change font heights
;res@tiXAxisFontHeightF  = 0.015                    ; change font heights
;res@tiYAxisFontHeightF  = 0.015                    ; change font heights

;res@gsnLeftString   = "dx = 50m; dz = 50m; dt = 0.1s"               ; add the gsn titles
res@gsnLeftString   = "dx = 100m; dz = 100m; dt = 0.1s"               ; add the gsn titles
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
	plots(0) = gsn_csm_contour(wks,data(i,:,:),res)
end do

;=================================
end
;=================================
