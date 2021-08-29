from tkinter import *
import json, re, numpy
import tkinter.ttk as ttk
import tkinter.messagebox as msgbox
from matplotlib import pyplot as plt
import matplotlib.cm as cm

# 아미노산 json 파일 읽기
with open("과학탐구대회/amino_acid.json", "r", encoding="UTF8") as f:
    amino_data = json.load(f)

# RGB 이미지 제작을 위한 배열 제작 함수
def make_rgb_img(sequenceX, sequenceY, img_size):
    
    # 비교대상1
    codon_1, i = [], 0
    for codon in sequenceX: # 염기서열의 코돈
        if codon == [0]: # 리스트 값이 0이면 for문 나옴
            break
        elif str(codon[0]).isdigit(): # 0이 아닌 숫자면 for문 나옴
            break
        for nuc in codon: # A T G C 를 각각 1 2 3 4 로 치환
            if nuc == "a":
                codon_1.append(1)
            elif nuc == "t":
                codon_1.append(2)
            elif nuc == "g":
                codon_1.append(3)
            elif nuc == "c":
                codon_1.append(4)
        sequenceX.remove(codon)
        sequenceX.insert(i, codon_1)
        i+=1
        codon_1 = []

    Img_x = [] # 이미지화할 배열
    codon_2, codon_3 = [], []
    for codon in sequenceX: # 염기서열의 코돈
        for i in range(3):
            if codon == [0]: # 리스트 값이 0이면 흰색으로 지정
                for j in range(3):
                    codon_2.append(255)
                break
            if codon[i] == 3: # 보고서에 작성한 공식 연산 과정
                N = 16*9*codon[int((i+1)%3)]*codon[int((i+2)%3)]/codon[int(i%3)]**2 -1
            else:
                N = 16*codon[int((i+1)%3)]*codon[int((i+2)%3)]/codon[int(i%3)]**2 -1
            codon_2.append(int(N))
        codon_3.append(codon_2)
        codon_2 = []
        if len(codon_3) == img_size:
            Img_x.append(codon_3) # 리스트에 추가
            codon_3= [] 

    # 비교대상2
    codon_1, i = [], 0
    for codon in sequenceY:
        if codon == [0]:
            break
        elif str(codon[0]).isdigit():
            break
        for nuc in codon:
            if nuc == "a":
                codon_1.append(1)
            elif nuc == "t":
                codon_1.append(2)
            elif nuc == "g":
                codon_1.append(3)
            elif nuc == "c":
                codon_1.append(4)
        sequenceY.remove(codon)
        sequenceY.insert(i, codon_1)
        i+=1
        codon_1 = []

    Img_y = [] 
    codon_2, codon_3 = [], []
    for codon in sequenceY:
        for i in range(3):
            if codon == [0]:
                for j in range(3):
                    codon_2.append(255)
                break
            if codon[i] == 3:
                N = 16*9*codon[int((i+1)%3)]*codon[int((i+2)%3)]/codon[int(i%3)]**2 -1
            else:
                N = 16*codon[int((i+1)%3)]*codon[int((i+2)%3)]/codon[int(i%3)]**2 -1
            codon_2.append(int(N))
        codon_3.append(codon_2)
        codon_2 = []
        if len(codon_3) == img_size:
            Img_y.append(codon_3)
            codon_3= []

    fig = plt.figure(figsize=(8,8))
    
    # 비교대상1 염기서열 시각화
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.imshow(Img_x)
    ax1.set_title('Target_X')
    ax1.axis("off")

    # 비교대상2 염기서열 시각화
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.imshow(Img_y)
    ax2.set_title('Target_Y')
    ax2.axis("off")

    plt.show()

# cmap 이미지 제작 함수
def make_cmap_img(sequenceX, sequenceY, img_size):
    num_x, num_y = [], []
    color_X, color_Y = [], []

    # 입력받은 염기서열 코돈에 대응되는 아미노산 찾아 색 지정
    for codon in sequenceX: # 염기서열의 코돈
        if codon == [0]: # 리스트 값이 0이면 색 지정 값도 0으로 지정
            num_x.append(0)        
        for i in range(21): # 아미노산 json 파일에서 같은 것을 찾아
            for j in range(len(amino_data['amino_acid'][i]['codon'])):
                if amino_data['amino_acid'][i]['codon'][j] == codon:
                    num_x.append(amino_data['amino_acid'][i]['N'] + j*10)
                if len(num_x) == img_size:
                    color_X.append(num_x)
                    num_x = [] # 각 코돈에 해당하는 값을 색 지정을 위한 리스트에 저장

    for codon in sequenceY:
        if codon == [0]:
            num_y.append(0)        
        for i in range(21):
            for j in range(len(amino_data['amino_acid'][i]['codon'])):
                if amino_data['amino_acid'][i]['codon'][j] == codon:
                    num_y.append(amino_data['amino_acid'][i]['N'] + j*10)
                if len(num_y) == img_size:
                    color_Y.append(num_y)
                    num_y = []

    # 이미지화 할 배열
    Img_x = numpy.array(color_X)
    Img_y = numpy.array(color_Y)

    # 이미지 색 계열 선택
    color = cmb_color.get()
    if color == "Grey":
        color = cm.Greys
    elif color == "Purple":
        color = cm.Purples
    elif color == "Blue":
        color = cm.Blues
    elif color == "Green":
        color = cm.Greens
    elif color == "Orange":
        color = cm.Oranges
    elif color == "Red":
        color = cm.Reds

    fig = plt.figure(figsize=(8,8))

    # 비교대상1 염기서열 시각화
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.imshow(Img_x, cmap=color)
    ax1.set_title('Target_X')
    ax1.axis("off")

    # 비교대상2 염기서열 시각화
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.imshow(Img_y, cmap=color)
    ax2.set_title('Target_Y')
    ax2.axis("off")

    plt.show()

# 이미지 크기 결정 함수
def find_size(sequenceX, sequenceY):
    len_seqX, len_seqY = len(sequenceX), len(sequenceY) # 염기서열 길이 변수에 저장

    # 비교할 두 리스트 중 더 긴 것의 길이를 변수에 저장
    if len_seqX > len_seqY: 
        find_size = len_seqX
    else:
        find_size = len_seqY

    # 제곱수 찾기
    for size in range(101): 
        if size**2 >= find_size :
            img_size=size
            break

    # 배열의 남은 공간을 0으로 채움
    for i in range(img_size**2):
        if len_seqX < img_size**2:
            sequenceX.append([0])
        if len_seqY < img_size**2:
            sequenceY.append([0])

    # RGB, Colormap 중 선택한 것 실행        
    if method_var.get() == "RGB": # RGB 이미지 제작 함수 실행
        make_rgb_img(sequenceX, sequenceY, img_size)
    elif method_var.get() == "cmap": # Colormap 이미지 제작 함수 실행
        make_cmap_img(sequenceX, sequenceY, img_size)

# 염기서열 정렬 함수
def sort():
    # 입력받은 염기서열에서 공백과 숫자 제거
    targetX = txt_targetX.get("1.0", END).replace(" ", "").replace("\n", "")
    targetX = re.sub(r'[0-9]+', '', targetX)
    targetY = txt_targetY.get("1.0", END).replace(" ", "").replace("\n", "")
    targetY = re.sub(r'[0-9]+', '', targetY)

    sequenceX, sequenceY = [], [] # 리스트로 정리된 염기서열

    # 코돈별로 구분해 이중  리스트로 정리
    for i in range(int(len(targetX)/3)):
        targetX_1 = list(targetX[i*3] + targetX[i*3+1] + targetX[i*3+2])
        sequenceX.append(targetX_1)

    for i in range(int(len(targetY)/3)):
        targetY_1 = list(targetY[i*3] + targetY[i*3+1] + targetY[i*3+2])
        sequenceY.append(targetY_1)

    find_size(sequenceX, sequenceY) # 이미지 크기 결정 함수 실행

# input창에 염기서열이 입력되었는지 확인
def start():
    if len(txt_targetX.get("1.0", END)) == 1 or len(txt_targetY.get("1.0", END)) == 1:
        msgbox.showwarning("경고", "염기서열을 입력하세요")
        return # 염기서열을 입력하지 않았으면 경고창 출력
    else:
        sort() # 염기서열을 입력했으면 염기서열 정렬 함수 실행

root = Tk()
root.title("염기서열 시각화")
root.geometry("800x380")
root.resizable(False, False)

# 염기서열 입력 프레임
frame_input = Frame(root, relief="solid", bd=3)
frame_input.pack(fill="both", padx=5, pady=5)

# 비교 대상 1
frame_targetX = LabelFrame(frame_input, text="비교대상1", width=60)
frame_targetX.pack(side="left", fill="both", padx=5, pady=5)

scrollbarX = Scrollbar(frame_targetX)
scrollbarX.pack(side="right", fill="y")

txt_targetX = Text(frame_targetX, width=50, height=15, yscrollcommand=scrollbarX.set)
txt_targetX.pack(padx=5, pady=5)
scrollbarX.config(command=txt_targetX.yview)

# 비교 대상 2
frame_targetY = LabelFrame(frame_input, text="비교대상2", width=60)
frame_targetY.pack(side="right", fill="both",  padx=5, pady=5)

scrollbarY = Scrollbar(frame_targetY)
scrollbarY.pack(side="right", fill="y")

txt_targetY = Text(frame_targetY, width=50, height=15, yscrollcommand=scrollbarY.set)
txt_targetY.pack(padx=5, pady=5)
scrollbarY.config(command=txt_targetY.yview)

# 옵션 프레임
frame_option = LabelFrame(root, text="옵션")
frame_option.pack(fill="x", padx=5, pady=5)

# 표현방식 선택 레이블
lbl_method = Label(frame_option, text="표현방식", width=8)
lbl_method.pack(side="left", padx=5, pady=5)

# 표현방식 선택 콤보
method_var = StringVar()
btn_rgb = Radiobutton(frame_option, text="RGB", value="RGB", variable=method_var)
btn_rgb.select()
btn_colormap = Radiobutton(frame_option, text="cmap", value="cmap", variable=method_var)
btn_rgb.pack(side="left")
btn_colormap.pack(side="left")

# cmap 색상 선택 콤보
opt_color = ["Grey", "Purple", "Blue", "Green", "Orange", "Red"]
cmb_color = ttk.Combobox(frame_option, state="readonly", value=opt_color, width=10)
cmb_color.current(0)
cmb_color.pack(side="right", padx=5, pady=5)

# cmap 색상 선택 레이블
lbl_color = Label(frame_option, text="Colormap을 골랐다면 색상을 선택하세요.")
lbl_color.pack(side="right", padx=5, pady=5)

# 실행 프레임
frame_run = Frame(root)
frame_run.pack(fill="x", padx=5, pady=5)

# 시작 버튼
btn_close = Button(frame_run, padx=5, pady=5, text="닫기", width=12, command=root.quit)
btn_close.pack(side="right", padx=5, pady=5)

# 종료 버튼
btn_start = Button(frame_run, padx=5, pady=5, text="시작", width=12, command=start)
btn_start.pack(side="right", padx=5, pady=5)

root.mainloop()
