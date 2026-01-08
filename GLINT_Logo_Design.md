# GLINT Logo Design Concepts

## 🎨 Design Philosophy

**Core Concept**: "闪光时刻" - The moment of discovery when molecular glue interfaces are revealed

**Visual Elements**:
- **Glint/Sparkle**: Represents insight, discovery, illumination
- **Ternary Structure**: E3-Glue-POI triangle
- **Interface**: Connection, interaction, binding

---

## 🥇 Recommended Design: "Ternary Sparkle"

### Visual Description
```
        E3
       /  \
      /    \
     /  ✦   \    ← Central sparkle (glue molecule)
    /        \
  POI ——————— Glue
```

### Design Specifications

**Icon Version (64x64 px)**:
- Three nodes forming triangle (E3, POI, Glue)
- Central sparkle/star symbol
- Clean lines, modern geometric style
- Color: Gradient blue (#3b82f6 → #60a5fa)

**Full Logo (Horizontal)**:
```
  ✦
 / \     GLINT
/   \    Glue Interface Analyzer
```

**Color Palette**:
- Primary: Blue (#3b82f6) - Trust, science, technology
- Accent: Light Blue (#60a5fa) - Innovation, clarity
- Sparkle: White/Yellow (#fbbf24) - Discovery, insight
- Text: Dark Gray (#1f2937) - Professional, readable

### File Formats Needed
1. **Icon**: 
   - `glint_icon.png` (64x64, 128x128, 256x256)
   - `glint_icon.svg` (vector, scalable)
   - `glint_icon.ico` (Windows)

2. **Logo**:
   - `glint_logo_horizontal.png` (300 DPI)
   - `glint_logo_horizontal.svg`
   - `glint_logo_vertical.png`

3. **Banner**:
   - `glint_banner.png` (1200x400, for GitHub)

---

## 🥈 Alternative Design 1: "Interface Wave"

### Visual Description
```
[Protein 1] ≈≈✦≈≈ [Protein 2]
```

**Concept**: Ripple effect from glue molecule creating interface

**Pros**: 
- Dynamic, shows interaction
- Emphasizes "interface" concept

**Cons**: 
- More complex to render at small sizes
- Less distinctive

---

## 🥉 Alternative Design 2: "G-Letter Fusion"

### Visual Description
```
  ╭─╮
 │ ● │  ← G letter + molecular structure
  ╰─╯
   |
  ✦
```

**Concept**: Letter G integrated with molecular representation

**Pros**: 
- Clear brand identity
- Memorable shape

**Cons**: 
- Less directly related to "interface" concept
- May look generic

---

## 📐 Technical Implementation

### SVG Code for Recommended Design (Simplified)

```svg
<svg width="64" height="64" xmlns="http://www.w3.org/2000/svg">
  <defs>
    <linearGradient id="blueGradient" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" style="stop-color:#3b82f6;stop-opacity:1" />
      <stop offset="100%" style="stop-color:#60a5fa;stop-opacity:1" />
    </linearGradient>
  </defs>
  
  <!-- Triangle (Ternary Complex) -->
  <line x1="32" y1="10" x2="10" y2="50" stroke="url(#blueGradient)" stroke-width="2"/>
  <line x1="32" y1="10" x2="54" y2="50" stroke="url(#blueGradient)" stroke-width="2"/>
  <line x1="10" y1="50" x2="54" y2="50" stroke="url(#blueGradient)" stroke-width="2"/>
  
  <!-- Nodes (E3, POI, Glue) -->
  <circle cx="32" cy="10" r="4" fill="#3b82f6"/>
  <circle cx="10" cy="50" r="4" fill="#3b82f6"/>
  <circle cx="54" cy="50" r="4" fill="#3b82f6"/>
  
  <!-- Central Sparkle (Glue Interface) -->
  <path d="M32,25 L34,30 L39,30 L35,33 L37,38 L32,35 L27,38 L29,33 L25,30 L30,30 Z" 
        fill="#fbbf24" stroke="#f59e0b" stroke-width="0.5"/>
</svg>
```

### PNG Generation (Python + Pillow)

```python
from PIL import Image, ImageDraw
import math

def create_glint_icon(size=256):
    img = Image.new('RGBA', (size, size), (255, 255, 255, 0))
    draw = ImageDraw.Draw(img)
    
    # Triangle vertices
    e3 = (size//2, size//5)
    poi = (size//5, size*4//5)
    glue = (size*4//5, size*4//5)
    
    # Draw triangle
    draw.line([e3, poi], fill=(59, 130, 246), width=size//32)
    draw.line([e3, glue], fill=(59, 130, 246), width=size//32)
    draw.line([poi, glue], fill=(59, 130, 246), width=size//32)
    
    # Draw nodes
    r = size//16
    draw.ellipse([e3[0]-r, e3[1]-r, e3[0]+r, e3[1]+r], fill=(59, 130, 246))
    draw.ellipse([poi[0]-r, poi[1]-r, poi[0]+r, poi[1]+r], fill=(59, 130, 246))
    draw.ellipse([glue[0]-r, glue[1]-r, glue[0]+r, glue[1]+r], fill=(59, 130, 246))
    
    # Draw sparkle (center)
    center = (size//2, size//2)
    draw_sparkle(draw, center, size//8, (251, 191, 36))
    
    return img

def draw_sparkle(draw, center, size, color):
    """Draw a 4-point star sparkle"""
    cx, cy = center
    points = []
    for i in range(8):
        angle = math.pi * i / 4
        r = size if i % 2 == 0 else size // 3
        x = cx + r * math.cos(angle)
        y = cy + r * math.sin(angle)
        points.append((x, y))
    draw.polygon(points, fill=color, outline=(245, 158, 11))

# Generate icons
for size in [64, 128, 256]:
    img = create_glint_icon(size)
    img.save(f'glint_icon_{size}.png')
```

---

## 🎯 Usage Guidelines

### Minimum Size
- Icon: 16x16 px (must remain recognizable)
- Logo with text: 120px width minimum

### Clear Space
- Maintain padding equal to 1/4 of icon height on all sides

### Color Variations
1. **Full Color**: Primary use (blue gradient + yellow sparkle)
2. **Monochrome**: Single blue (#3b82f6) for print
3. **White**: For dark backgrounds
4. **Black**: For light backgrounds, formal documents

### Don'ts
- ❌ Don't distort aspect ratio
- ❌ Don't change colors arbitrarily
- ❌ Don't add effects (shadows, glows) except sparkle
- ❌ Don't place on busy backgrounds

---

## 📦 Deliverables Checklist

- [ ] `glint_icon_64.png`
- [ ] `glint_icon_128.png`
- [ ] `glint_icon_256.png`
- [ ] `glint_icon.svg`
- [ ] `glint_icon.ico`
- [ ] `glint_logo_horizontal.png`
- [ ] `glint_logo_horizontal.svg`
- [ ] `glint_logo_vertical.png`
- [ ] `glint_banner_github.png` (1200x400)
- [ ] `glint_favicon.ico` (16x16, 32x32)

---

## 🚀 Next Steps

1. **Generate icon files** using provided Python script
2. **Create SVG master file** for scalability
3. **Update PyMOL plugin** to use new icon
4. **Update GitHub repository** with new branding
5. **Create favicon** for documentation site

Would you like me to generate the actual icon files now?

