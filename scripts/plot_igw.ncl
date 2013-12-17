load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;=================================
begin
;=================================
; Settings.
;---------------------------------
nx = 600
nz = 200

;nx = 300
;nz = 100

;nx = 150
;nz = 50

;nx = 120
;nz = 40

;nt = 1
;nt = 3
;nt = 4
;nt = 10
;nt = 11
;nt = 16
;nt = 31
;nt = 65
;nt = 83
;nt = 101
;nt = 188
;nt = 161
;nt = 201
nt = 301
;nt = 311
;nt = 401
;nt = 1001
var = "theta_1"
;var = "theta"
;var = "pi_1"
;var = "u"
;var = "w"
halo = 4
;=================================
;Read data.
;---------------------------------
setfileoption ("bin", "ReadByteOrder", "BigEndian")
;data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz,nx/),"float")
data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz+2*halo,nx+2*halo/),"double")
data@_FillValue = -9999.
;print(data)
;last = data(nt-1,:,:)
;=================================
; Plot.
;---------------------------------
wks = gsn_open_wks ("pdf","igw")                  ; open workstation
;gsn_define_colormap(wks,"temp_diff_18lev")
;gsn_define_colormap(wks,"WhiteYellowOrangeRed")
;gsn_define_colormap(wks,"WhiteBlue")
gsn_define_colormap(wks,"BlueWhiteOrangeRed")

res = True

res@gsnMaximize = True
;res@gsnPaperOrientation = "portrait"
res@vpWidthF  = 1.0                    ; change the aspect ratio
res@vpHeightF = 0.3
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
;res@gsnSpreadColorStart = 8
;res@gsnSpreadColorEnd = 60

res@cnLevelSelectionMode = "ManualLevels"   ; set explicit contour levels
;res@cnMinLevelValF          = -0.0015
;res@cnMaxLevelValF          = 0.003

res@cnMinLevelValF          = -0.0030
res@cnMaxLevelValF          = 0.0030
res@cnLevelSpacingF         = 0.0005

;res@cnFillColors = (/7,8,9,10,12,13,14,15,16,17,18,19,20/)
;res@cnFillColors = (/8,9,10,12,13,14,15,16,17,18,19,20/)
;res@cnFillColors = (/5,6,7,8,9,10,12,13,14,15,16,17,18,19,20/)

; Line
;res@cnFillOn = False
;res@cnLinesOn = True
;res@cnLineLabelsOn = True
;res@cnMonoLevelFlag = False
;res@cnLevelFlags = "LineAndLabel"

res@tiMainString = "Inertia Gravity Waves"
;res@tiXAxisString      = "X"  
;res@tiYAxisString      = "Z (height)"

res@tiMainFontHeightF  = 0.015                    ; change font heights
;res@tiXAxisFontHeightF  = 0.015                    ; change font heights
;res@tiYAxisFontHeightF  = 0.015                    ; change font heights

;res@gsnLeftString   = "dx = 500m; dz = 50m; dt = 0.1s"               ; add the gsn titles
res@gsnLeftString   = "dx = 1000m; dz = 100m; dt = 0.1s"               ; add the gsn titles
;res@gsnLeftString   = "dx = 2500m; dz = 250m; dt = 0.1s"               ; add the gsn titles
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
;tt = 299
;last = data(tt,:,:)
;res@gsnRightStringFontHeightF   = 0.015
;res@gsnRightString  = "Step: " + tt*100
;plots(0) = gsn_csm_contour(wks,last,res)

;=================================
end
;=================================
