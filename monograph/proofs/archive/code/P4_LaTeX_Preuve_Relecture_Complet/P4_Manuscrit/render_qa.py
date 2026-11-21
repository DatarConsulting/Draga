from pathlib import Path
import fitz,json,math
from PIL import Image,ImageDraw
root=Path(__file__).resolve().parent
doc=fitz.open(root/'P4_Manuscrit_Complet.pdf')
qa=root/'qa';qa.mkdir(exist_ok=True)
out=[];bounds=[]
for i,page in enumerate(doc):
 pix=page.get_pixmap(matrix=fitz.Matrix(.65,.65),alpha=False)
 img=Image.frombytes('RGB',[pix.width,pix.height],pix.samples)
 img.thumbnail((290,412))
 out.append(img)
 for block in page.get_text('dict')['blocks']:
  if block['type']!=0:continue
  for line in block['lines']:
   for sp in line['spans']:
    x0,y0,x1,y1=sp['bbox']
    if x0<10 or y0<7 or x1>page.rect.width-10 or y1>page.rect.height-7:
     bounds.append({'page':i+1,'text':sp['text'][:100],'bbox':sp['bbox']})
for start in range(0,len(out),9):
 sheet=Image.new('RGB',(930,1340),'#dce3e5');draw=ImageDraw.Draw(sheet)
 for j,img in enumerate(out[start:start+9]):
  x=15+(j%3)*305;y=22+(j//3)*443
  draw.text((x,y-16),f'PDF page {start+j+1}',fill='black')
  sheet.paste(img,(x,y))
 sheet.save(qa/f'contact_{start//9+1:02d}.jpg',quality=90)
(qa/'layout_bounds.json').write_text(json.dumps({'pages':len(doc),'out_of_page_text':bounds},indent=2))
print({'pages':len(doc),'contact_sheets':math.ceil(len(doc)/9),'out_of_page_text':len(bounds)})
