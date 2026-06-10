from fpdf import FPDF
from PIL import Image
import os

def create_combined_pdf(plot_dir="../../plots", output_pdf="../../results/combined_plots.pdf"):
    if not os.path.exists(plot_dir):
        print(f"Directory {plot_dir} does not exist.")
        return

    pdf = FPDF()

    images = []
    for root, dirs, files in os.walk(plot_dir):
        for f in files:
            if f.lower().endswith(('.png', '.jpg', '.jpeg')):
                images.append(os.path.join(root, f))

    images.sort() # Ensure consistent order

    if not images:
        print(f"No images found in {plot_dir}")
        return

    print(f"Found {len(images)} images. creating PDF...")

    for image_path in images:
        try:
            pdf.add_page()

            rel_path = os.path.relpath(image_path, plot_dir)

            with Image.open(image_path) as img:
                w, h = img.size

                max_w = 190
                max_h = 277

                aspect = h / w

                if aspect > (max_h / max_w):
                     print_h = max_h
                     print_w = max_h / aspect
                else:
                     print_w = max_w
                     print_h = max_w * aspect

                pdf.image(image_path, x=10, y=10, w=print_w, h=print_h)
                pdf.set_font('Arial', '', 10)
                pdf.text(10, 290, f"{rel_path}") # Print filename at bottom
        except Exception as e:
            print(f"Could not process {image_path}: {e}")

    pdf.output(output_pdf)
    print(f"Successfully created {output_pdf}")

if __name__ == "__main__":
    create_combined_pdf()
