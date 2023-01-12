import numpy as np
from PIL import Image
import math
import shapely.geometry.polygon as shplyg
import random

"""---------超参数设置-----------"""
SEQUENT = 1
RANDOM = 0
X0, Y0 = 0, 0  # 原点坐标
mult = 0.7  # 初始化边框乘的比例系数
width_step = 30  # 扩展边框宽度的步长
height_step = 20  # 扩展边框高度的步长
a = 2  # 重叠长度                  # 适应度函数比例系数
b = 1  # 角重叠占比
c = 2.5  # 拼接后图片面积/包围图片的最小矩形面积
d = 1  # 包围图片的最小矩形的长宽比
"""-------------------------------"""


# 矩形类的定义
class Rect():
    def __init__(self, x, y, height, width, id=None):
        self.x = x
        self.y = y  # 左上角
        self.height = height  # 高
        self.width = width  # 宽
        self.id = id  # 编号

        self.area = self.width * self.height  # 面积
        self.center = (self.x + self.width / 2, self.y + self.height / 2)  # 中心点坐标
        self.diaglen = np.sqrt(self.height ** 2 + self.width ** 2)  # 对角线长度
        self.pntlist = [(self.x, self.y), (self.x, self.y + self.height),  # 从左上角顶点开始逆时针
                        (self.x + self.width, self.y + self.height),
                        (self.x + self.width, self.y)]  # 顶点列表
        pntLst = []
        for i in range(4):
            pntLst.append(self.pntlist[i])
        self.pntTuple = tuple(pntLst)

        # 更新顶点信息

    def updatepnts(self):
        self.pntlist = [(self.x, self.y), (self.x, self.y + self.height),
                        (self.x + self.width, self.y + self.height), (self.x + self.width, self.y)]

        # 更新矩形的位置

    def updatePosition(self, x, y):
        self.x = x
        self.y = y
        self.updatepnts()

        # 取左上角顶点和右下角顶点

    def getCorner(self):
        lu = self.pntlist[0]
        rd = self.pntlist[2]
        return lu, rd


# 多边形类的定义
class Polygon():
    def __init__(self, pntlist):
        self.pntlist = pntlist  # 顶点坐标
        self.nospnt = len(self.pntlist)
        self.lu = (0, 0)  # 拼接成的多边形默认顶点为左上角顶点，即原点
        self.rd = np.max(self.pntlist, 0)  # 把多边形围起来的矩阵的右上角坐标
        self.center = np.mean([self.lu, self.rd], 0)  # 多边形的中心点


# 计算包围当前已拼接图片的最小矩形的长和宽
def minBox(Rect_list):
    LUx, LUy = X0, Y0
    RDx, RDy = -10 ** 6, -10 ** 6
    for rect in Rect_list:
        rd = rect.getCorner()[1]
        rdx, rdy = rd
        if rdx > RDx:
            RDx = rdx
        if rdy > RDy:
            RDy = rdy
    height = RDy - LUy
    width = RDx - LUx
    return [width, height]


# 初始化 将图片转化为对应的矩形
def init_image(image_list):
    """
    输入图片列表， 返回矩形列表
    """
    Rect_list = []
    area_list = []
    height, width, count, org_x, org_y, flag = 0, 0, 0, 0, 0, 0
    for file in image_list:
        image = Image.open(file)
        height = image.size[1]  # 图片高度&矩形高度
        width = image.size[0]  # 图片宽度&矩形宽度
        area_list.append(height * width)
        if flag:
            org_x += width
            flag = 0
        else:
            org_y += height
            flag = 1
        Rect_list.append(Rect(org_x, org_y, id=count, height=height, width=width))
        count += 1
    return Rect_list, area_list


# 对图片/矩形按面积从大到小排序
def sortbyArea(Rect_list, area_list):
    """
    返回排序好的矩形列表
    """
    sorted_Rect_list = []
    index_list = sorted(range(len(area_list)), key=lambda i: area_list[i], reverse=True)
    for i in range(len(Rect_list)):
        sorted_Rect_list.append(Rect_list[index_list[i]])
    return sorted_Rect_list


# 更新角点
def updatepntlst(key_point, keep_pntlst, move_pntlst):  # key_point为连接的关键点
    # print('---------------------------------------------------',keep_pntlst,move_pntlst,key_point)
    keep_key_index = keep_pntlst.index(key_point)
    move_key_index = move_pntlst.index(key_point)
    # keep
    if (keep_key_index == 0):  # 相邻角点
        last_keep_point = keep_pntlst[len(keep_pntlst) - 1]
    else:
        last_keep_point = keep_pntlst[keep_key_index - 1]
    if (keep_key_index == (len(keep_pntlst) - 1)):
        next_keep_point = keep_pntlst[0]
    else:
        next_keep_point = keep_pntlst[keep_key_index + 1]
    # move
    if (move_key_index == 0):  # 相邻角点
        last_move_point = move_pntlst[len(move_pntlst) - 1]
    else:
        last_move_point = move_pntlst[move_key_index - 1]
    if (move_key_index == (len(move_pntlst) - 1)):
        next_move_point = move_pntlst[0]
    else:
        next_move_point = move_pntlst[move_key_index + 1]

    keep_pntlst.pop(keep_key_index)
    move_pntlst.pop(move_key_index)

    slice1 = keep_pntlst[0:keep_pntlst.index(last_keep_point) + 1]
    slice2 = keep_pntlst[keep_pntlst.index(next_keep_point):]
    if (keep_pntlst.index(next_keep_point) == 0):
        slice2 = []
    total_pntlst = []
    total_pntlst = slice1 + move_pntlst + slice2
    if (last_move_point == next_keep_point):
        total_pntlst.pop(total_pntlst.index(next_keep_point))
        total_pntlst.pop(total_pntlst.index(next_keep_point))
    if (next_move_point == last_keep_point):
        total_pntlst.pop(total_pntlst.index(last_keep_point))
        total_pntlst.pop(total_pntlst.index(last_keep_point))
    # print(slice1, slice2, '\n返回值', total_pntlst)
    return total_pntlst


# 更新拼接之后已拼接图片&多边形的关键点
# def updateKeypnt(CombinedImg, MOVE):
#     # 输入已拼接好图片(Polygon类)和当前拼接的图片，输出更新后的Polygon
#     # comb_point = MOVE.pntlist[0]
#     temp_MOVE_pntlist = []
#     temp_Combined_pntlist = []
#     same_pntlist = []
#     for point in MOVE.pntlist:
#         if point in CombinedImg.pntlist:
#             same_pntlist.append(point)
#     for point in CombinedImg.pntlist:
#         if point not in same_pntlist:
#             temp_Combined_pntlist.append(point)
#     for point in MOVE.pntlist:
#         if point not in same_pntlist:
#             temp_MOVE_pntlist.append(point)
#     temp_point_list = temp_Combined_pntlist + temp_MOVE_pntlist
#     sorted_point_list = [temp_point_list.pop(0)]
#     # print(temp_point_list)
#     # print(sorted_point_list)
#
#     flag = 0
#     while len(temp_point_list):
#         point_same_line = []
#         for point in temp_point_list:
#             if point[flag] == sorted_point_list[-1][flag]:
#                 point_same_line.append(point)
#         sorted_point_list.append(tuple(np.min(point_same_line, 0)))
#         temp_point_list.remove(tuple(np.min(point_same_line, 0)))
#         if flag == 0:
#             flag += 1
#         else:
#             flag -= 1
#     CombinedImg = Polygon(sorted_point_list)
#     return CombinedImg, sorted_point_list


# 根据所有图片总面积计算初始化边框大小   边框定义为矩形类！

def computeFrame(area_list, mult):
    total_area = np.sum(area_list)
    print(f'The total area of the images is {total_area}.')
    width = math.ceil(mult * np.sqrt(total_area))
    frame_rect = Rect(X0, Y0, width=width, height=width)
    return frame_rect


# 检测当前的边框长宽比是否处于0.9-1.1
def checkWHratio(frame_rect):
    if 0.9 <= frame_rect.width / frame_rect.height <= 1.1:
        print(f'The width height ratio now is {frame_rect.width / frame_rect.height}.')
        return True
    else:
        print(f'The width height ratio now is {frame_rect.width / frame_rect.height}!!!')
        return False


# 扩展当前的边框
def extendFrame(frame_rect):
    new_frame_rect = Rect(frame_rect.x, frame_rect.y, frame_rect.width + width_step,
                          frame_rect.height + height_step)
    if checkWHratio(new_frame_rect):
        return new_frame_rect
    else:
        new_frame_rect = Rect(frame_rect.x, frame_rect.y, frame_rect.width + width_step,
                              (frame_rect.width + width_step) * 0.9)
        return new_frame_rect


# 输出当前状态信息
def print_info(Combing_Img_list, Rect_list, MOVE):
    print(f'Now {len(Combing_Img_list)} images have been combined, there are {len(Rect_list) + 1} '
          f'images need to combined\n Process...{MOVE.id} {image_list[MOVE.id]}.')
    # print('-'*50)


# 计算当前拼接结果的适应度
def computeFitness(MOVE, Combined_Img):
    poly1 = shplyg.Polygon(Combined_Img.pntlist)
    poly2 = shplyg.Polygon(MOVE.pntlist)
    # 计算重叠长度
    total_length = poly1.length + poly2.length
    common_area = poly1.intersection(poly2)  # 区域是一条线
    length = common_area.length
    length_ratio = length / total_length
    # print(length_ratio)

    # 计算角重叠占比
    num = 0.0
    for i in range(len(MOVE.pntlist)):
        move_point = MOVE.pntlist[i]
        for j in range(len(Combined_Img.pntlist)):
            keep_point = Combined_Img.pntlist[j]
            if (move_point == keep_point):
                num += 1
    total_num = len(Combined_Img.pntlist) + len(MOVE.pntlist)
    point_ratio = num / total_num
    # print('----------', point_ratio)

    # 计算矩形实体占比
    body_area = poly1.area + poly2.area
    # total_pntlst = []
    total_pntlst = Combined_Img.pntlist + MOVE.pntlist
    ld = np.min(total_pntlst, 0)  # 把多边形围起来的矩阵的左下角坐标
    ru = np.max(total_pntlst, 0)  # 把多边形围起来的矩阵的右上角坐标
    rect_area = abs((ru[0] - ld[0]) * (ld[1] - ru[1]))
    area_ratio = body_area / rect_area
    # print('面积占比', area_ratio)

    # 拼接后包围已拼接图片的矩形长宽比，最后数值越大越好
    egde_ratio_0 = abs((ru[0] - ld[0]) / (ld[1] - ru[1]))
    egde_ratio = 1 / (abs(egde_ratio_0 - 1) + 1)
    # print(egde_ratio)

    fitness = a * length_ratio + b * point_ratio + c * area_ratio + d * egde_ratio
    return fitness


def isOverlap(shplgi, shplgj):
    itstn = shplgi.intersection(shplgj)
    # print(itstn.area, itstn.length)
    if itstn.area != 0:
        return True
    else:
        return False


# 检测当前拼接的MOVE图片&矩形是否和已经拼接好的图片重叠
def checkOverlap(MOVE, Combined_Img):
    check1 = shplyg.Polygon(Combined_Img.pntlist)  # 类型转换
    # if not check1.is_valid:
    #     print('检测失效！！！-----------')
    #     return INVALID
    check2 = shplyg.Polygon(MOVE.pntlist)
    return isOverlap(check1, check2)


# 检测当前拼接的MOVE图片&矩形是否超出边框
def checkOverframe(MOVE, frame_rect):
    lu, rd = MOVE.getCorner()
    lux, luy = lu[0], lu[1]
    rdx, rdy = rd[0], rd[1]
    if lux < frame_rect.x or luy < frame_rect.y \
            or rdx > frame_rect.x + frame_rect.width or rdy > frame_rect.x + frame_rect.height:
        # print('overFrame')
        return True
    else:
        # print('OK')
        return False


# 依次选取图片进行拼接
def Select(MOVE, Combined_Img, frame_rect):
    find, maxFitness = 0, 0
    ret_MOVE = MOVE
    for point in Combined_Img.pntlist:
        temp_MOVE = Rect(point[0], point[1], width=MOVE.width, height=MOVE.height, id=MOVE.id)
        if checkOverlap(temp_MOVE, Combined_Img):
            continue
        if checkOverframe(temp_MOVE, frame_rect):
            continue
        fitness = computeFitness(temp_MOVE, Combined_Img)
        find = 1
        if fitness > maxFitness:
            maxFitness = fitness
            ret_MOVE = temp_MOVE
    return find, maxFitness, ret_MOVE


# 输出拼接后的图片
def plotIMG(img_width, img_height, Combined_Img_list, image_list, output):
    out = Image.new('RGBA', (img_width, img_height))
    for rect in Combined_Img_list:
        # print('Process...', image_list[rect.id])
        image = Image.open(image_list[rect.id])
        out.paste(image, (rect.x, rect.y))
    out.save(output)


# 拼接图片函数
def ImgCombine(Rect_list, area_list, output):
    """---initial---"""
    # 初始化第一张图片的位置
    Rect_list[0].updatePosition(X0, Y0)
    # 取出第一张图片
    Combined_Img_list = []
    Combined_Img_list.append(Rect_list.pop(0))
    Combined_Img = Polygon(Combined_Img_list[0].pntlist)

    # 初始化边框
    frame_rect = computeFrame(area_list=area_list, mult=mult)
    print(f'--The initial frame_rect: area--{frame_rect.area}, '
          f'width--{frame_rect.width}, height--{frame_rect.height}.--')

    # 待拼接的图片&矩形
    """---依次拼接---"""
    while len(Rect_list):
        MOVE = Rect_list.pop(0)
        print_info(Combined_Img_list, Rect_list, MOVE)

        # 选取适应值最高的拼接状态进行拼接
        find, maxFitness, temp_MOVE = Select(MOVE, Combined_Img, frame_rect)

        # 所有状态都不满足条件，扩展边框
        if not find:
            """扩展边框"""
            # Rect_list.append(MOVE)
            Rect_list.insert(0, MOVE)
            frame_rect = extendFrame(frame_rect)
            print('The frame rect has been extended!!!')
            print(f'--The new frame_rect: area--{frame_rect.area}, '
                  f'width--{frame_rect.width}, height--{frame_rect.height}.--')
        else:
            print(f'The max fitness is {maxFitness}.')
            print('--' * 30)
            # 将MOVE图片加入已拼接的图片
            MOVE = temp_MOVE
            # Combined_Img, cmbd_keypnt = updateKeypnt(Combined_Img, MOVE)
            cmbd_keypntlst = updatepntlst((MOVE.x, MOVE.y), Combined_Img.pntlist, MOVE.pntlist)
            Combined_Img = Polygon(cmbd_keypntlst)
            Combined_Img_list.append(MOVE)

    # print(Combined_Img.pntlist)
    print('All images have been combined.')

    # 最终包围所有拼接好的图片的最小矩形的长宽
    minBox_width, minBox_height = Combined_Img.rd
    # print('The min box width-', minBox_width, '--The min box height- ', minBox_height)
    # print(f'The width height ratio is {minBox_width / minBox_height}.')

    # 判断最小矩形的宽高比是否在0.9-1.1
    if 0.9 <= minBox_width / minBox_height <= 1.1:
        img_width, img_height = minBox_width, minBox_height
        print('True')
    else:
        print('False')
        if minBox_width > minBox_height:
            img_width = minBox_width
            img_height = math.ceil(minBox_width / 1.1)
        if minBox_width < minBox_height:
            img_height = minBox_height
            img_width = math.ceil(0.9 * minBox_height)

    print(f'The image width is {img_width}, The image height is {img_height}.')
    print(f'The image width height ratio is {img_width / img_height}.')
    print(f'The image area is {img_width * img_height}.')
    # 输出图片
    plotIMG(img_width, img_height, Combined_Img_list, image_list, output)

    area = img_width * img_height
    return area


if __name__ == "__main__":
    # 对应的case
    case = 4
    output = 'case' + str(case) + '.png'
    filename = 'image_name_case' + str(case) + '.txt'

    # image_list 为储存每个case中图片文件文件名的列表
    image_list = []
    with open(filename) as file_obj:
        for line in file_obj:
            image_list.append(line.rstrip())

    # Rect_list 为对于的矩形列表
    Rect_list, area_list = init_image(image_list)
    print(image_list)
    print('--' * 10, 'Case', case, '--' * 10)
    print(f'There are {len(image_list)} need to combine.---')
    print('--' * 30)

    Mode = SEQUENT
    if Mode == SEQUENT:
        print('SEQUENT')
        Rect_list = sortbyArea(Rect_list, area_list)
    if Mode == RANDOM:
        print('RANDON')
        random.shuffle(Rect_list)

    print('--' * 30)
    print('Init........')
    print(f'Init frame mult:{mult}, extend frame width step: {width_step}, ',
          f'height_step:{height_step} ',
          f'fitness rate: a={a}, b={b}, c={c}, d={d}')

    Combined_Img_area = ImgCombine(Rect_list=Rect_list, area_list=area_list, output=output)
