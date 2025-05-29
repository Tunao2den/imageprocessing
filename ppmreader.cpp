#include<stdio.h>
#include<conio.h>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>
#include<string.h>
#include <fcntl.h>
#include <malloc.h>
#include <math.h>
#define PI 3.1415926535897932384626433832795
#pragma pack(1)

#pragma region structs

struct ppm_header
{
	char pgmtype1;
	char pgmtype2;
	int pwidth;
	int pheight;
	int pmax;
};
struct ppm_file
{
	struct ppm_header *pheader;
	unsigned char *rdata,*gdata,*bdata;
};

struct hsi_pixel {
    float h;  // Hue [0, 360]
    float s;  // Saturation [0, 1]
    float i;  // Intensity [0, 1]
};

struct hsi_image {
    int width;
    int height;
    struct hsi_pixel *data;  // Boyut: width * height
};

struct ycbcr_pixel {
    float y;   // Luma (Y)
    float cb;  // Blue-difference chroma (Cb)
    float cr;  // Red-difference chroma (Cr)
};

struct ycbcr_image {
    int width;
    int height;
    struct ycbcr_pixel *data;  // Boyut: width * height
};


#pragma endregion

void get_image_data(const char *filename,struct ppm_file *image);
void write_image(const char *filename,struct ppm_file *image);
void write_image(const char *filename,struct ppm_file *image)
{
	FILE *fp;
	int i,max=0;
	fp=fopen(filename,"wb");
	fputc(image->pheader->pgmtype1,fp);
	fputc(image->pheader->pgmtype2,fp);
	fputc('\n',fp);
	fprintf(fp,"%d %d\n",image->pheader->pwidth,image->pheader->pheight);
	fprintf(fp,"%d\n",255/*max*/);
	for(i=0;i<image->pheader->pwidth*image->pheader->pheight;i++)
	{
		fwrite(&image->rdata[i],1,1,fp);
		fwrite(&image->gdata[i],1,1,fp);
		fwrite(&image->bdata[i],1,1,fp);
	}
	fclose(fp);
}
void get_image_data(const char *filename, struct ppm_file *image ) 
{
	FILE* fp;
	int i=0;
	char temp[256];
	image->pheader=(struct ppm_header *)malloc(sizeof(struct ppm_header));
	fp = fopen(filename, "rb" );
	if (fp==NULL) 
	{
		printf("File is not opened: %s.\n\n", filename);
		exit(1);
	}
	printf ("The PPM File : %s...\n", filename);
	fscanf (fp, "%s", temp);
	if (strcmp(temp, "P6") == 0) 
	{
		image->pheader->pgmtype1=temp[0];
		image->pheader->pgmtype2=temp[1];
		fscanf (fp, "%s", temp);
		if (temp[0]=='#') 
		{
			while(fgetc(fp)!='\n');
			fscanf (fp, "%d %d\n",&image->pheader->pwidth,&image->pheader->pheight);
			fscanf (fp, "%d\n", &image->pheader->pmax);
		}
		else
		{
			sscanf (temp, "%d", &image->pheader->pwidth);
			fscanf (fp, "%d", &image->pheader->pheight);
			fscanf (fp, "%d", &image->pheader->pmax);
		}
		image->rdata=(unsigned char *)malloc(image->pheader->pheight*image->pheader->pwidth*sizeof(unsigned char));
		image->gdata=(unsigned char *)malloc(image->pheader->pheight*image->pheader->pwidth*sizeof(unsigned char));
		image->bdata=(unsigned char *)malloc(image->pheader->pheight*image->pheader->pwidth*sizeof(unsigned char));
		if (image->rdata==NULL) printf("Memory problem\n");
		for(i=0;i<image->pheader->pwidth*image->pheader->pheight;i++)
		{
			fread(&image->rdata[i],1,1,fp);
			fread(&image->gdata[i],1,1,fp);
			fread(&image->bdata[i],1,1,fp);
		}
	}
	else 
	{
		printf ("\nError! The file is not a PPM file");
		exit(1);
	}
	fclose(fp);
}

void copy_ppm_file(struct ppm_file *src, struct ppm_file *dest) {
    dest->pheader = (struct ppm_header *)malloc(sizeof(struct ppm_header));
    if (!dest->pheader) {
        printf("Memory allocation failed for header.\n");
        exit(1);
    }
    memcpy(dest->pheader, src->pheader, sizeof(struct ppm_header));

    int total_pixels = src->pheader->pwidth * src->pheader->pheight;

    dest->rdata = (unsigned char *)malloc(total_pixels * sizeof(unsigned char));
    dest->gdata = (unsigned char *)malloc(total_pixels * sizeof(unsigned char));
    dest->bdata = (unsigned char *)malloc(total_pixels * sizeof(unsigned char));

    if (!dest->rdata || !dest->gdata || !dest->bdata) {
        printf("Memory allocation failed for pixel data.\n");
        exit(1);
    }

    memcpy(dest->rdata, src->rdata, total_pixels);
    memcpy(dest->gdata, src->gdata, total_pixels);
    memcpy(dest->bdata, src->bdata, total_pixels);
}

struct hsi_image* ppm_to_hsi_image(struct ppm_file *ppm) {
    int width = ppm->pheader->pwidth;
    int height = ppm->pheader->pheight;
    int total_pixels = width * height;

    struct hsi_image *hsi_img = (struct hsi_image *)malloc(sizeof(struct hsi_image));
    hsi_img->width = width;
    hsi_img->height = height;
    hsi_img->data = (struct hsi_pixel *)malloc(sizeof(struct hsi_pixel) * total_pixels);

    for (int i = 0; i < total_pixels; i++) {
        float r = ppm->rdata[i] / 255.0f;
        float g = ppm->gdata[i] / 255.0f;
        float b = ppm->bdata[i] / 255.0f;

        float num = 0.5f * ((r - g) + (r - b));
        float den = sqrt((r - g) * (r - g) + (r - b) * (g - b));
        float theta = 0.0f;

        if (den != 0)
            theta = acosf(fminf(fmaxf(num / den, -1.0f), 1.0f)); // clamp to [-1,1]
        
        float H = 0.0f;
        if (b <= g)
            H = theta * (180.0f / M_PI); // radian to degree
        else
            H = (2 * M_PI - theta) * (180.0f / M_PI);

        float min_rgb = fminf(r, fminf(g, b));
        float sum = r + g + b;
        float S = (sum == 0) ? 0 : 1 - (3 * min_rgb / sum);
        float I = sum / 3;

        hsi_img->data[i].h = H; // 0 - 360
        hsi_img->data[i].s = S; // 0 - 1
        hsi_img->data[i].i = I; // 0 - 1
    }

    return hsi_img;
}

void histogram_equalize_hsi(struct hsi_image *hsi_img) {
    int total_pixels = hsi_img->width * hsi_img->height;
    int hist[256] = {0};
    int cdf[256] = {0};
    int i;

    for (i = 0; i < total_pixels; i++) {
        int intensity = (int)(hsi_img->data[i].i * 255.0f);
        if (intensity > 255) intensity = 255;
        if (intensity < 0) intensity = 0;
        hist[intensity]++;
    }

    cdf[0] = hist[0];
    for (i = 1; i < 256; i++) {
        cdf[i] = cdf[i - 1] + hist[i];
    }

    float cdf_min = 0;
    for (i = 0; i < 256; i++) {
        if (cdf[i] != 0) {
            cdf_min = cdf[i];
            break;
        }
    }

    float map[256];
    for (i = 0; i < 256; i++) {
        map[i] = ((float)(cdf[i] - cdf_min) / (total_pixels - cdf_min)) * 255.0f;
        if (map[i] < 0) map[i] = 0;
        if (map[i] > 255) map[i] = 255;
    }

    for (i = 0; i < total_pixels; i++) {
        int old_intensity = (int)(hsi_img->data[i].i * 255.0f);
        if (old_intensity > 255) old_intensity = 255;
        if (old_intensity < 0) old_intensity = 0;

        float new_intensity = map[old_intensity] / 255.0f;
        hsi_img->data[i].i = new_intensity;
    }
}

struct ppm_file *hsi_to_ppm_image(struct hsi_image *hsi_img) {
    int width = hsi_img->width;
    int height = hsi_img->height;
    int total_pixels = width * height;

    struct ppm_file *image = (struct ppm_file *)malloc(sizeof(struct ppm_file));
    image->pheader = (struct ppm_header *)malloc(sizeof(struct ppm_header));
    image->pheader->pgmtype1 = 'P';
    image->pheader->pgmtype2 = '6';
    image->pheader->pwidth = width;
    image->pheader->pheight = height;
    image->pheader->pmax = 255;

    image->rdata = (unsigned char *)malloc(total_pixels);
    image->gdata = (unsigned char *)malloc(total_pixels);
    image->bdata = (unsigned char *)malloc(total_pixels);

    for (int i = 0; i < total_pixels; i++) {
        float H = hsi_img->data[i].h;
        float S = hsi_img->data[i].s;
        float I = hsi_img->data[i].i;

        float R, G, B;

        if (S == 0.0f) {
            R = G = B = I;
        } else {
            float h_deg = H * 360.0f;
            float h_rad = h_deg * M_PI / 180.0f;

            if (h_deg >= 0 && h_deg < 120) {
                B = I * (1 - S);
                R = I * (1 + (S * cos(h_rad)) / cos(M_PI / 3 - h_rad));
                G = 3 * I - (R + B);
            } else if (h_deg >= 120 && h_deg < 240) {
                h_rad = (h_deg - 120.0f) * M_PI / 180.0f;
                R = I * (1 - S);
                G = I * (1 + (S * cos(h_rad)) / cos(M_PI / 3 - h_rad));
                B = 3 * I - (R + G);
            } else {
                h_rad = (h_deg - 240.0f) * M_PI / 180.0f;
                G = I * (1 - S);
                B = I * (1 + (S * cos(h_rad)) / cos(M_PI / 3 - h_rad));
                R = 3 * I - (G + B);
            }
        }

        int r = (int)(R * 255.0f);
        int g = (int)(G * 255.0f);
        int b = (int)(B * 255.0f);

        if (r < 0) r = 0; if (r > 255) r = 255;
        if (g < 0) g = 0; if (g > 255) g = 255;
        if (b < 0) b = 0; if (b > 255) b = 255;

        image->rdata[i] = (unsigned char)r;
        image->gdata[i] = (unsigned char)g;
        image->bdata[i] = (unsigned char)b;
    }

    return image;
}


struct ycbcr_image *ppm_to_ycbcr_image(struct ppm_file *image) {
    int width = image->pheader->pwidth;
    int height = image->pheader->pheight;
    int total_pixels = width * height;

    struct ycbcr_image *ycbcr_img = (struct ycbcr_image *)malloc(sizeof(struct ycbcr_image));
    ycbcr_img->width = width;
    ycbcr_img->height = height;
    ycbcr_img->data = (struct ycbcr_pixel *)malloc(sizeof(struct ycbcr_pixel) * total_pixels);

    for (int i = 0; i < total_pixels; i++) {
        unsigned char r = image->rdata[i];
        unsigned char g = image->gdata[i];
        unsigned char b = image->bdata[i];

        float y  = 0.299f * r + 0.587f * g + 0.114f * b;
        float cb = -0.168736f * r - 0.331264f * g + 0.5f * b + 128.0f;
        float cr = 0.5f * r - 0.418688f * g - 0.081312f * b + 128.0f;

        ycbcr_img->data[i].y = y;
        ycbcr_img->data[i].cb = cb;
        ycbcr_img->data[i].cr = cr;
    }

    return ycbcr_img;
}



void histogram_equalize_ycbcr(struct ycbcr_image *img) {
    int total_pixels = img->width * img->height;
    int histogram[256] = {0};
    int cdf[256] = {0};

    for (int i = 0; i < total_pixels; i++) {
        int y_val = (int)(img->data[i].y); 
        if (y_val < 0) y_val = 0;
        if (y_val > 255) y_val = 255;
        histogram[y_val]++;
    }

    cdf[0] = histogram[0];
    for (int i = 1; i < 256; i++) {
        cdf[i] = cdf[i - 1] + histogram[i];
    }

    int cdf_min = 0;
    for (int i = 0; i < 256; i++) {
        if (cdf[i] != 0) {
            cdf_min = cdf[i];
            break;
        }
    }

    for (int i = 0; i < total_pixels; i++) {
        int y_val = (int)(img->data[i].y);
        if (y_val < 0) y_val = 0;
        if (y_val > 255) y_val = 255;

        int y_new = round(((float)(cdf[y_val] - cdf_min) / (total_pixels - cdf_min)) * 255.0f);
        img->data[i].y = (float)y_new;
    }
}


void YCbCr_to_rgb(struct ycbcr_image *ycbcr_img, struct ppm_file *ppm_img) {
    int width = ycbcr_img->width;
    int height = ycbcr_img->height;
    int total_pixels = width * height;

    ppm_img->pheader = (struct ppm_header *)malloc(sizeof(struct ppm_header));
    ppm_img->pheader->pgmtype1 = 'P';
    ppm_img->pheader->pgmtype2 = '6';
    ppm_img->pheader->pwidth = width;
    ppm_img->pheader->pheight = height;
    ppm_img->pheader->pmax = 255;

    ppm_img->rdata = (unsigned char *)malloc(total_pixels);
    ppm_img->gdata = (unsigned char *)malloc(total_pixels);
    ppm_img->bdata = (unsigned char *)malloc(total_pixels);

    for (int i = 0; i < total_pixels; i++) {
        float Y  = ycbcr_img->data[i].y;
        float Cb = ycbcr_img->data[i].cb;
        float Cr = ycbcr_img->data[i].cr;

        int R = round(Y + 1.402 * (Cr - 128));
        int G = round(Y - 0.344136 * (Cb - 128) - 0.714136 * (Cr - 128));
        int B = round(Y + 1.772 * (Cb - 128));

        if (R < 0) R = 0; if (R > 255) R = 255;
        if (G < 0) G = 0; if (G > 255) G = 255;
        if (B < 0) B = 0; if (B > 255) B = 255;

        ppm_img->rdata[i] = (unsigned char)R;
        ppm_img->gdata[i] = (unsigned char)G;
        ppm_img->bdata[i] = (unsigned char)B;
    }
}


double calculate_snr(struct ppm_file *original, struct ppm_file *processed) {
    int total_pixels = original->pheader->pwidth * original->pheader->pheight;
    double signal = 0.0;
    double noise = 0.0;

    for (int i = 0; i < total_pixels; i++) {
        // Red
        double o_r = (double)original->rdata[i];
        double p_r = (double)processed->rdata[i];
        signal += o_r * o_r;
        noise += (o_r - p_r) * (o_r - p_r);

        // Green
        double o_g = (double)original->gdata[i];
        double p_g = (double)processed->gdata[i];
        signal += o_g * o_g;
        noise += (o_g - p_g) * (o_g - p_g);

        // Blue
        double o_b = (double)original->bdata[i];
        double p_b = (double)processed->bdata[i];
        signal += o_b * o_b;
        noise += (o_b - p_b) * (o_b - p_b);
    }

    if (noise == 0) {
        return INFINITY; // perfect match
    }

    return 10.0 * log10(signal / noise);
}

double calculate_psnr(struct ppm_file *original, struct ppm_file *processed) {
    int total_pixels = original->pheader->pwidth * original->pheader->pheight;
    double mse = 0.0;

    for (int i = 0; i < total_pixels; i++) {
        // Red
        double diff_r = (double)original->rdata[i] - (double)processed->rdata[i];
        mse += diff_r * diff_r;

        // Green
        double diff_g = (double)original->gdata[i] - (double)processed->gdata[i];
        mse += diff_g * diff_g;

        // Blue
        double diff_b = (double)original->bdata[i] - (double)processed->bdata[i];
        mse += diff_b * diff_b;
    }

    mse /= (3.0 * total_pixels);

    if (mse == 0.0) {
        return INFINITY; // perfect match
    }

    double max_pixel = 255.0;
    return 10.0 * log10((max_pixel * max_pixel) / mse);
}


void save_image_data(const char *filename, struct ppm_file *image) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        printf("Could not open file for writing: %s\n", filename);
        exit(1);
    }

    fprintf(fp, "%c%c\n", image->pheader->pgmtype1, image->pheader->pgmtype2);
    fprintf(fp, "%d %d\n", image->pheader->pwidth, image->pheader->pheight);
    fprintf(fp, "%d\n", image->pheader->pmax);

    int total_pixels = image->pheader->pwidth * image->pheader->pheight;

    for (int i = 0; i < total_pixels; i++) {
        fwrite(&image->rdata[i], 1, 1, fp);
        fwrite(&image->gdata[i], 1, 1, fp);
        fwrite(&image->bdata[i], 1, 1, fp);
    }

    fclose(fp);
}

#include <stdio.h>
#include <stdlib.h>

void save_ppm_as_bmp(struct ppm_file *image, const char *filename) {
    FILE *f;
    int width = image->pheader->pwidth;
    int height = image->pheader->pheight;
    int padding = (4 - (width * 3) % 4) % 4;
    int filesize = 54 + (3 * width + padding) * height;

    unsigned char bmpfileheader[14] = {
        'B','M',  // BMP header
        0,0,0,0,  // Dosya boyutu (4 byte)
        0,0,      // reserved
        0,0,      // reserved
        54,0,0,0  // header size (14 + 40)
    };

    unsigned char bmpinfoheader[40] = {
        40,0,0,0,      // info header size
        0,0,0,0,       // width
        0,0,0,0,       // height
        1,0,           // planes
        24,0,          // bits per pixel
        0,0,0,0,       // compression
        0,0,0,0,       // image size (can be 0)
        0x13,0x0B,0,0, // X pixels per meter (2835)
        0x13,0x0B,0,0, // Y pixels per meter (2835)
        0,0,0,0,       // total colors
        0,0,0,0        // important colors
    };

    bmpfileheader[2] = (unsigned char)(filesize    );
    bmpfileheader[3] = (unsigned char)(filesize>> 8);
    bmpfileheader[4] = (unsigned char)(filesize>>16);
    bmpfileheader[5] = (unsigned char)(filesize>>24);

    bmpinfoheader[4] = (unsigned char)(width    );
    bmpinfoheader[5] = (unsigned char)(width>> 8);
    bmpinfoheader[6] = (unsigned char)(width>>16);
    bmpinfoheader[7] = (unsigned char)(width>>24);

    bmpinfoheader[8] = (unsigned char)(height    );
    bmpinfoheader[9] = (unsigned char)(height>> 8);
    bmpinfoheader[10]= (unsigned char)(height>>16);
    bmpinfoheader[11]= (unsigned char)(height>>24);

    f = fopen(filename, "wb");
    fwrite(bmpfileheader, 1, 14, f);
    fwrite(bmpinfoheader, 1, 40, f);

    for (int y = height - 1; y >= 0; y--) {
        for (int x = 0; x < width; x++) {
            int i = y * width + x;
            unsigned char bgr[3] = {
                image->bdata[i],
                image->gdata[i],
                image->rdata[i]
            };
            fwrite(bgr, 1, 3, f);
        }
        // padding
        for (int p = 0; p < padding; p++)
            fputc(0x00, f);
    }

    fclose(f);
    //printf("Saved as BMP: %s\n", filename);
}

void open_image(const char *filename) {
    char command[256];
	snprintf(command, sizeof(command), "start %s", filename);  // Windows
    system(command);
}

void free_ppm_file(struct ppm_file *image) {
    if (image->pheader) free(image->pheader);
    if (image->rdata) free(image->rdata);
    if (image->gdata) free(image->gdata);
    if (image->bdata) free(image->bdata);
}

void free_hsi_image(struct hsi_image *image) {
    if (image) {
        if (image->data) {
            free(image->data);
        }
        free(image);
    }
}

void free_ycbcr_image(struct ycbcr_image *image) {
    if (image) {
        if (image->data) {
            free(image->data);
        }
        free(image);
    }
}

int main()
{	
	struct ppm_file picture;
	int i;
	get_image_data("mandrill.ppm",&picture);

	struct ppm_file picture1, picture2;
	copy_ppm_file(&picture, &picture1);
	copy_ppm_file(&picture, &picture2);

	printf("pgmtype...=%c%c\n",picture.pheader->pgmtype1,picture.pheader->pgmtype2);
	printf("width...=%d\n",picture.pheader->pwidth);
	printf("height...=%d\n",picture.pheader->pheight);
	printf("max gray level...=%d\n",picture.pheader->pmax);

	struct hsi_image *hsi_pic1 = ppm_to_hsi_image(&picture1);
	histogram_equalize_hsi(hsi_pic1);
	struct ppm_file *picture1_eq = hsi_to_ppm_image(hsi_pic1);


	struct ycbcr_image *ycbcr_pic2 = ppm_to_ycbcr_image(&picture2);
	histogram_equalize_ycbcr(ycbcr_pic2);
	YCbCr_to_rgb(ycbcr_pic2, &picture2);

	double snr_hsi = calculate_snr(&picture, picture1_eq);
	double snr_ycbcr = calculate_snr(&picture, &picture2);

	printf("SNR (HSI Equalized): %.2f dB\n", snr_hsi);
	printf("SNR (YCbCr Equalized): %.2f dB\n", snr_ycbcr);


	double psnr_hsi = calculate_psnr(&picture, picture1_eq);
	double psnr_ycbcr = calculate_psnr(&picture, &picture2);

	printf("PSNR (HSI Equalized): %.2f dB\n", psnr_hsi);
	printf("PSNR (YCbCr Equalized): %.2f dB\n", psnr_ycbcr);


	save_image_data("original.ppm", &picture);
	save_image_data("equalized_hsi.ppm", picture1_eq);
	save_image_data("equalized_ycbcr.ppm", &picture2);

	save_ppm_as_bmp(&picture, "original.bmp");
	save_ppm_as_bmp(picture1_eq, "hsi_equalized.bmp");
	save_ppm_as_bmp(&picture2, "ycbcr_equalized.bmp");


	open_image("original.bmp");
	open_image("hsi_equalized.bmp");
	open_image("ycbcr_equalized.bmp");

	write_image("pnr.ppm",&picture);

	free_ppm_file(&picture);
	free_ppm_file(picture1_eq);
	free_ppm_file(&picture2);
	free_hsi_image(hsi_pic1);
	free_ycbcr_image(ycbcr_pic2);

	return 0;
} 