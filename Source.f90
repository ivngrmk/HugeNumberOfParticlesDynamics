program PlotBar
    use IFQWIN
    ! ! !Всё для графики! ! !
    integer(4), parameter :: uTERMINAL = 7
    integer(4), parameter :: uPLOT = 8
    type(wxycoord)     :: PlotXY
    type(qwinfo)      :: FrameSize
    type(qwinfo)      :: TerminalSize
    type(qwinfo)      :: PlotSize
    type(windowconfig):: TerminalWindow
    type(windowconfig):: PlotWindow
    integer(2)        :: i2Status
    integer(4)        :: i4Status
    integer(4)        :: i4LetterX, i4LetterY
    integer(4)        :: i4Columns, i4Rows
    logical(4) l4Status
    logical(2) l2Status
    ! ! !Частицы! ! !
    integer, parameter :: NumberOfParticles = 12
    type particle
        real x
        real y
        real vx
        real vy
    end type particle
    type vector
        real x
        real y
    end type
    type(particle) oldParticle(12)
    type(particle) newParticle(12)
    real, parameter :: vmax = 1.0
    real impulsex, impulsey
    ! ! !Среда! ! !
    real, parameter :: Lx = 8
    real, parameter :: Ly = 8
    ! ! !Счётчики! ! !
    integer i, j, k
    ! ! !Параметры симуляции! ! !
    real, parameter :: dt = 0.02
    ! ! !Используется в динамике! ! !
    type(vector) partAccelerationOld
    type(vector) partAccelerationNew
    !InitPlotWindow(uPLOTinternal)
    real xl, yl, xr, yr, scale_width, x_scale, y_scale
    real x, y
    integer(4) uPLOTinternal
    !NewAcceleration(partNumb) и OldAcceleration(partNumb)
    type(vector) r
    real r_length
    !InitParticles()
    real gridStepx, gridStepy
    integer gridPartNumbx, gridPartNumby
    real randomNumber
    !StepBehaviour()
    logical(2) l2Exit
    integer(4) mouse_x, mouse_y
    integer(4) keystate
    ! ! ! ! ! ! ! ! ! ! ! ! !
    call RANDOM_SEED()
    call InitWindow()
    !allocate(oldParticle(NumberOfParicles))
    !allocate(newParticle(NumberOfParicles))
    call InitParticles()
    call VisualizeState()
    call StepBehaviour()
    contains
    subroutine StepBehaviour()
        l2Exit = .TRUE.
        i4Status = SETACTIVEQQ(uTERMINAL)
        i4Status = WAITONMOUSEEVENT(MOUSE$LBUTTONDOWN,keystate,mouse_x,mouse_yy)
        do while (l2Exit == .TRUE.)
            if (keystate == MOUSE$KS_RBUTTON) then
                l2Exit = .FALSE.
            else
                do k=1,1
                    call SleepQQ(10)
                    call SimulationStep()
                end do
                call VisualizeState()
            end if
        end do
    end subroutine StepBehaviour
    ! ! ! ! ! ! ! ! ! ! ! ! !
    subroutine VisualizeState()
        real radius
        i4Status = SETACTIVEQQ(uPLOT)
        l2Status = InitPlotWindow(uPLOT)
        l2Status = SetColor(4)
        radius = min(Lx,Ly) / 30.0
        do i =1,NumberOfParticles
            i2Status = ELLIPSE_W($GFILLINTERIOR, oldParticle(i).x - radius, oldParticle(i).y + radius, oldParticle(i).x + radius, oldParticle(i).y - radius)
        end do
    end subroutine VisualizeState
    ! ! ! ! ! ! ! ! ! ! ! ! !
    subroutine InitParticles()
        gridPartNumbx = 3;
        if (DBLE(NumberOfParticles) / gridPartNumbx /= INT(DBLE(NumberOfParticles) / gridPartNumbx) ) then
            write(uTERMINAL,*) "Error while making grid of particles."
        else
            gridPartNumby = NumberOfParticles / gridPartNumbx
        end if
        gridStepx = floor(Lx / 2.0 / gridPartNumbx); gridStepy = floor(Ly / 2.0 / gridPartNumby)
        do i=1,gridPartNumbx
            do j=1,gridPartNumby
                oldParticle( (i-1)*gridPartNumby + j ).x = i * gridStepx
                oldParticle( (i-1)*gridPartNumby + j ).y = j * gridStepy
                call RANDOM_NUMBER(randomNumber)
                oldParticle( (i-1)*gridPartNumby + j ).vx = vmax*(2*randomNumber-1)
                call RANDOM_NUMBER(randomNumber)
                randomNumber = randomNumber - 0.5
                if (randomNumber < 0) then
                    randomNumber = -1.0
                else
                    randomNumber = 1.0
                end if
                oldParticle( (i-1)*gridPartNumby + j ).vy = randomNumber * sqrt(vmax**2 - oldParticle( (i-1)*gridPartNumby + j ).vx**2)
            end do
        end do
        impulsex = 0.0; impulsey = 0.0
        do i=1,NumberOfParticles - 1
            impulsex = impulsex + oldParticle(i).vx
            impulsey = impulsey + oldParticle(i).vy
        end do
        oldParticle(NumberOfParticles).vx = -impulsex
        oldParticle(NumberOfParticles).vy = -impulsey
    end subroutine InitParticles
    ! ! ! ! ! ! ! ! ! ! ! ! !
    subroutine SimulationStep()
        do i=1,NumberOfParticles
            partAccelerationOld = OldAcceleration(i)
            newParticle(i).x = oldParticle(i).x + oldParticle(i).vx * dt + 1.0 / 2.0 * partAccelerationOld.x*dt**2
            newParticle(i).y = oldParticle(i).y + oldParticle(i).vy * dt + 1.0 / 2.0 * partAccelerationOld.y*dt**2
        end do
        do i=1,NumberOfParticles
            partAccelerationOld = OldAcceleration(i)
            partAccelerationNew = NewAcceleration(i)
            newParticle(i).vx = oldParticle(i).vx + 1.0 / 2.0 * (partAccelerationOld.x + partAccelerationNew.x) * dt
            newParticle(i).vy = oldParticle(i).vy + 1.0 / 2.0 * (partAccelerationOld.y + partAccelerationNew.y) * dt
        end do
        do i=1,NumberOfParticles
            oldParticle(i).x = newParticle(i).x
            oldParticle(i).y = newParticle(i).y
            if (oldParticle(i).x > Lx) then
                oldParticle(i).x = newParticle(i).x - Lx
            end if
            if (oldParticle(i).x < 0) then
                oldParticle(i).x = newParticle(i).x + Lx
            end if
            if (oldParticle(i).y > Ly) then
                oldParticle(i).y = newParticle(i).y - Ly
            end if
            if (oldParticle(i).y < 0) then
                oldParticle(i).y = newParticle(i).y + Ly
            end if
            oldParticle(i).vx = newParticle(i).vx
            oldParticle(i).vy = newParticle(i).vy
        end do
    end subroutine SimulationStep
    ! ! ! ! ! ! ! ! ! ! ! ! !
    type(vector) function OldAcceleration(partNumb)
        integer partNumb
        OldAcceleration.x = 0.0; OldAcceleration.y = 0.0;
        do j = 1, NumberOfParticles
            if (j /= partNumb) then
                r.x = oldParticle(j).x - oldParticle(partNumb).x
                r.y = oldParticle(j).y - oldParticle(partNumb).y
                if (oldParticle(j).x - oldParticle(partNumb).x > Lx/2.0) then
                    r.x = oldParticle(j).x - oldParticle(partNumb).x - Lx
                end if
                if (oldParticle(j).x - oldParticle(partNumb).x < -Lx/2.0) then
                    r.x = oldParticle(j).x - oldParticle(partNumb).x + Lx
                end if
                if (oldParticle(j).y - oldParticle(partNumb).y > Ly/2.0) then
                    r.y = oldParticle(j).y - oldParticle(partNumb).y - Ly
                end if
                if (oldParticle(j).y - oldParticle(partNumb).y < -Ly/2.0) then
                    r.y = oldParticle(j).y - oldParticle(partNumb).y + Ly
                end if
                r_length = sqrt(r.x**2 + r.y**2)
                OldAcceleration.x = OldAcceleration.x + 4*(6 / r_length**8 - 12 / r_length**14)*r.x
                OldAcceleration.y = OldAcceleration.y + 4*(6 / r_length**8 - 12 / r_length**14)*r.y
            end if
        end do
    end function OldAcceleration
    ! ! ! ! ! ! ! ! ! ! ! ! !
    type(vector) function NewAcceleration(partNumb)
        integer partNumb
        NewAcceleration.x = 0.0; NewAcceleration.y = 0.0;
        do j = 1, NumberOfParticles
            if (j /= partNumb) then
                r.x = newParticle(j).x - newParticle(partNumb).x
                r.y = newParticle(j).y - newParticle(partNumb).y
                if (newParticle(j).x - newParticle(partNumb).x > Lx/2.0) then
                    r.x = newParticle(j).x - newParticle(partNumb).x - Lx
                end if
                if (newParticle(j).x - newParticle(partNumb).x < -Lx/2.0) then
                    r.x = newParticle(j).x - newParticle(partNumb).x + Lx
                end if
                if (newParticle(j).y - newParticle(partNumb).y > Ly/2.0) then
                    r.y = newParticle(j).y - newParticle(partNumb).y - Ly
                end if
                if (newParticle(j).y - newParticle(partNumb).y < -Ly/2.0) then
                    r.y = newParticle(j).y - newParticle(partNumb).y + Ly
                end if
                r_length = sqrt(r.x**2 + r.y**2)
                NewAcceleration.x = NewAcceleration.x + 4*(6 / r_length**8 - 12 / r_length**14)*r.x
                NewAcceleration.y = NewAcceleration.y + 4*(6 / r_length**8 - 12 / r_length**14)*r.y
            end if
        end do
    end function NewAcceleration
    ! ! ! ! ! ! ! ! ! ! ! ! !
    logical(2) function InitPlotWindow(uPLOTinternal) !Строит в выбранном окне декартову систему координат.
        integer uPLOTinternal
        l2Status = SETACTIVEQQ(uPLOTinternal)
        call SetViewPort(0, PlotWindow.numypixels, PlotWindow.numxpixels, 0)
        l2Status = SetBkColor(15); call ClearScreen($GVIEWPORT); !Окраска всего экрана в белый.
        xl = -1.0; yl = -1.0; xr = 9.0; yr = 9.0; scale_width = 0.1; x_scale = 1.0; y_scale = 1.0 !Обязательно должны содержать начало координат.
        l2Status = SetWindow(.TRUE., DBLE(xl), DBLE(yl), DBLE(xr), DBLE(yr))
        l2Status = SetColor(0)
        x = xl
        do while (ceiling(x) <= floor(xr)) !Градуировка шкалы абсцисс.
            call MoveTo_w(DBLE(ceiling(x)), DBLE(0.0 - scale_width), PlotXY)
            l2Status = LineTo_w(DBLE(ceiling(x)), DBLE(0.0 + scale_width))
            x = x + x_scale
        end do
        y = yl
        do while (ceiling(y) <= floor(yr)) !Градуировка шкалы ординат.
            call MoveTo_w(DBLE(0.0 - scale_width), DBLE(ceiling(y)), PlotXY)
            l2Status = LineTo_w(DBLE(0.0 + scale_width), DBLE(ceiling(y)))
            y = y + y_scale
        end do
        l2Status = SetColor(4) !Рисуем сами оси.
        call MoveTo_w(DBLE(xl), DBLE(0.0), PlotXY)
        l2Status = LineTo_w(DBLE(xr), DBLE(0.0))
        call MoveTo_w(DBLE(0.0), DBLE(yl), PlotXY)
        l2Status = LineTo_w(DBLE(0.0), DBLE(yr))
    end function InitPlotWindow
    ! ! ! ! ! ! ! ! ! ! ! ! !
subroutine InitWindow()
    !Размер шрифта.
    i4LetterX=8
    i4LetterY=16
    !Установка параметров внешнего окна.
    !Размер.
    i2Status=GETWSIZEQQ(QWIN$FRAMEWINDOW,QWIN$SIZEMAX,FrameSize)
    FrameSize.H=3*FrameSize.H/4 !! set window height 3/4 of max height
    FrameSize.W=3*FrameSize.W/4 !! set window width 3/4 of max width
    i2Status=SETWSIZEQQ(QWIN$FRAMEWINDOW,FrameSize)
    i4Columns = FrameSize.W / i4LetterX - 1 !Выражение ширины внешнего окна в единицах символов.
    i4Rows = FrameSize.H / i4LetterY        !Выражение высоты внешнего окна в единицах символов.
    !Внутренние окна.
    !Окно под консоль.
    open(uTERMINAL,file='user',title="Terminal")
    i4Status = SETACTIVEQQ(uTERMINAL)
    i4Status = GETWINDOWCONFIG(TerminalWindow)
    i2Status=GETWSIZEQQ(uTERMINAL,QWIN$SIZEMAX,TerminalSize)
    !Окно под график.
    open(uPLOT,file='user',title="Plot")
    i4Status = SETACTIVEQQ(uPLOT)
    i4Status = GETWINDOWCONFIG(PlotWindow)
    i2Status=GETWSIZEQQ(uPLOT,QWIN$SIZEMAX,PlotSize)
    !Настройка этих окон.
    TerminalSize.W = i4Columns / 2 - i4Columns / 16; PlotSize.W = i4Columns / 2 - i4Columns / 16;
    TerminalSize.H = i4Rows - i4Rows / 8;            PlotSize.H = i4Rows - i4Rows / 8;
    TerminalSize.X = 0;                              PlotSize.X = i4Columns / 2;
    i2Status=SETWSIZEQQ(uTERMINAL,TerminalSize)
    i2Status=SETWSIZEQQ(uPLOT,PlotSize)
    PlotWindow.numxpixels = PlotSize.W * i4LetterX
    PlotWindow.numypixels = PlotSize.H * i4LetterY
    i4Status = SETWINDOWCONFIG(PlotWindow)
    write(uTERMINAL,*) "Windows were successfully launched!"
    ! ! ! ! ! ! ! ! ! ! ! ! !
    l2Status = InitPlotWindow(uPLOT)
    write(uTERMINAL,*) "Plot was successfully launched!"
end subroutine
    ! ! ! ! ! ! ! ! ! ! ! ! !
end