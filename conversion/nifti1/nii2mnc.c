#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <limits.h>
#include <float.h>
#include <minc.h>
#include <ParseArgv.h>
#include <volume_io.h>

#include <time_stamp.h>
#include "nifti1_io.h"

/* Names of MINC spatial dimensions, in our "standard" world ordering.
 */
static const char *spatial_names[VIO_N_DIMENSIONS] = {
    MIxspace, MIyspace, MIzspace
};

void test_xform(mat44 m, int i, int j, int k)
{
    double x, y, z;

    x = m.m[VIO_X][0] * i + m.m[VIO_X][1] * j + m.m[VIO_X][2] * k
      + m.m[VIO_X][3];
    y = m.m[VIO_Y][0] * i + m.m[VIO_Y][1] * j + m.m[VIO_Y][2] * k
      + m.m[VIO_Y][3];
    z = m.m[VIO_Z][0] * i + m.m[VIO_Z][1] * j + m.m[VIO_Z][2] * k
      + m.m[VIO_Z][3];

    printf("%d %d %d => ", i, j, k);
    printf("%f %f %f\n", x, y, z);
}

static int usage(void)
{
    static const char msg[] = {
        "nii2mnc: Convert NIfTI-1 files to MINC format\n"
        "usage: nii2mnc [options] filename.nii [filename.mnc]\n"
    };
    fprintf(stderr, "%s", msg);
    return (-1);
}

static void
find_data_range(int datatype, size_t nvox, void *data, double range[2])
{
    size_t i;

    range[0] = DBL_MAX;
    range[1] = -DBL_MAX;

    for (i = 0; i < nvox; i++) {
        double tmp;

        switch (datatype) {
        case DT_INT8:
            tmp = (double) ((char *)data)[i];
            break;
        case DT_UINT8:
            tmp = (double) ((unsigned char *)data)[i];
            break;
        case DT_INT16:
            tmp = (double) ((short *)data)[i];
            break;
        case DT_UINT16:
            tmp = (double) ((unsigned short *)data)[i];
            break;
        case DT_INT32:
            tmp = (double) ((int *)data)[i];
            break;
        case DT_UINT32:
            tmp = (double) ((unsigned int *)data)[i];
            break;
        case DT_FLOAT32:
            tmp = (double) ((float *)data)[i];
            break;
        case DT_FLOAT64:
            tmp = (double) ((double *)data)[i];
            break;
        default:
            fprintf(stderr, "Data type %d not handled\n", datatype);
            break;
        }
        if (tmp < range[0]) {
            range[0] = tmp;
        }
        if (tmp > range[1]) {
            range[1] = tmp;
        }
    }
}

int
main(int argc, char **argv)
{
    /* NIFTI stuff */
    nifti_image *nii_ptr;

    /* MINC stuff */
    int mnc_fd;                 /* MINC file descriptor */
    nc_type img_mtype;          /* MINC memory data type */
    int img_msign;              /* MINC !0 if signed data */
    static nc_type img_vtype;   /* MINC voxel data type */
    static int img_vsign;       /* MINC !0 if signed data */
    int img_ndims;              /* MINC image dimension count */
    int img_dimids[MAX_VAR_DIMS]; /* MINC image dimension identifiers */
    int mnc_iid;                /* MINC Image variable ID */
    long img_start[MAX_VAR_DIMS]; /* MINC data starts */
    long img_count[MAX_VAR_DIMS]; /* MINC data counts */
    double img_vrange[2];       /* MINC valid min/max */
    double img_rrange[2];       /* MINC image min/max */
    double time_step;
    double time_start;
    int trivial_axes[VIO_N_DIMENSIONS] = {VIO_X, VIO_Y, VIO_Z};
    int spatial_axes[VIO_N_DIMENSIONS] = {VIO_X, VIO_Y, VIO_Z};
    double dim_starts[VIO_N_DIMENSIONS];
    double dim_steps[VIO_N_DIMENSIONS];
    double dim_dircos[VIO_N_DIMENSIONS][VIO_N_DIMENSIONS];
    VIO_Transform mnc_xfm;
    VIO_General_transform mnc_linear_xfm;
    mat44 nii_xfm;

    /* Other stuff */
    char out_str[1024];         /* Big string for filename */
    int i;                      /* Generic loop counter the first */
    int j;                      /* Generic loop counter the second */
    char *str_ptr;              /* Generic ASCIZ string pointer */
    int r;                      /* Result code. */
    static int qflag = 0;       /* Quiet flag (default is non-quiet) */
    static int rflag = 1;       /* Scan range flag */
    static char *mnc_ordered_dim_names[VIO_N_DIMENSIONS];

    static ArgvInfo argTable[] = {
        {"-byte", ARGV_CONSTANT, (char *) NC_BYTE, (char *)&img_vtype,
         "Write voxel data in 8-bit integer format."},
        {"-short", ARGV_CONSTANT, (char *) NC_SHORT, (char *)&img_vtype,
         "Write voxel data in 16-bit integer format."},
        {"-int", ARGV_CONSTANT, (char *) NC_INT, (char *)&img_vtype,
         "Write voxel data in 32-bit integer format."},
        {"-float", ARGV_CONSTANT, (char *) NC_FLOAT, (char *)&img_vtype,
         "Write voxel data in 32-bit floating point format."},
        {"-double", ARGV_CONSTANT, (char *) NC_DOUBLE, (char *)&img_vtype,
         "Write voxel data in 64-bit floating point format."},
        {"-signed", ARGV_CONSTANT, (char *) 1, (char *)&img_vsign,
         "Write integer voxel data in signed format."},
        {"-unsigned", ARGV_CONSTANT, (char *) 0, (char *)&img_vsign,
         "Write integer voxel data in unsigned format."},
        {"-noscanrange", ARGV_CONSTANT, (char *) 0, (char *)&rflag,
         "Do not scan data range."},
        {"-quiet", ARGV_CONSTANT, (char *) 0, (char *)&qflag,
         "Quiet operation"},
        {NULL, ARGV_END, NULL, NULL, NULL}
    };


    set_ncopts(0);      /* Clear global netCDF error reporting flag */

    img_vtype = NC_NAT;

    for (i = 0; i < MAX_VAR_DIMS; i++) {
      img_start[i] = 0;
    }

    if (ParseArgv(&argc, argv, argTable, 0)) {
        return usage();
    }
    if (argc < 2) {
        fprintf(stderr, "Too few arguments\n");
        return usage();
    }
    else if (argc == 2) {
        strcpy(out_str, argv[1]);
        str_ptr = strrchr(out_str, '.');
        if (str_ptr != NULL) {
            if (!strcmp(str_ptr, ".nii") || !strcmp(str_ptr, ".hdr")) {
                *str_ptr = '\0';
                strcat(out_str, ".mnc");
            }
        }
    }
    else if (argc == 3) {
        strcpy(out_str, argv[2]);
    }
    else {
        fprintf(stderr, "Extra arguments provided.\n");
        return usage();
    }

    /* Read in the entire NIfTI file. */
    nii_ptr = nifti_image_read(argv[1], 1);

    if (!qflag) {
        nifti_image_infodump(nii_ptr);
    }

    /* Determine the data type for output. */

    switch (nii_ptr->datatype) {
    case DT_INT8:
        img_msign = 1;
        img_mtype = NC_BYTE;
        img_vrange[0] = CHAR_MIN;
        img_vrange[1] = CHAR_MAX;
        break;
    case DT_UINT8:
        img_msign = 0;
        img_mtype = NC_BYTE;
        img_vrange[0] = 0;
        img_vrange[1] = UCHAR_MAX;
        break;
    case DT_INT16:
        img_msign = 1;
        img_mtype = NC_SHORT;
        img_vrange[0] = SHRT_MIN;
        img_vrange[1] = SHRT_MAX;
        break;
    case DT_UINT16:
        img_msign = 0;
        img_mtype = NC_SHORT;
        img_vrange[0] = 0;
        img_vrange[1] = USHRT_MAX;
        break;
    case DT_INT32:
        img_msign = 1;
        img_mtype = NC_INT;
        img_vrange[0] = INT_MIN;
        img_vrange[1] = INT_MAX;
        break;
    case DT_UINT32:
        img_msign = 0;
        img_mtype = NC_INT;
        img_vrange[0] = 0;
        img_vrange[1] = UINT_MAX;
        break;
    case DT_FLOAT32:
        img_msign = 1;
        img_mtype = NC_FLOAT;
        img_vrange[0] = -FLT_MAX;
        img_vrange[1] = FLT_MAX;
        break;
    case DT_FLOAT64:
        img_msign = 1;
        img_mtype = NC_DOUBLE;
        img_vrange[0] = -DBL_MAX;
        img_vrange[1] = DBL_MAX;
        break;
    default:
        fprintf(stderr, "Data type %d not handled\n", nii_ptr->datatype);
        break;
    }

    if (img_vtype == NC_NAT) {
        img_vsign = img_msign;
        img_vtype = img_mtype;
    }

    /* Calculate the starts, steps, and direction cosines. This only
     * be done properly if the file is NIfTI-1 file.  If it is an Analyze
     * file we have to resort to other methods...
     */
    if (nii_ptr->nifti_type != 0 &&
        (nii_ptr->sform_code != NIFTI_XFORM_UNKNOWN ||
         nii_ptr->qform_code != NIFTI_XFORM_UNKNOWN)) {

      if (nii_ptr->sform_code != NIFTI_XFORM_UNKNOWN) {
          if (!qflag) {
              printf("Using s-form transform:\n");
          }
          for (i = 0; i < 4; i++) {
              for (j = 0; j < 4; j++) {
                  nii_xfm.m[i][j] = nii_ptr->sto_xyz.m[i][j];
                  if (!qflag) {
                      printf("%8.4f, ", nii_ptr->sto_xyz.m[i][j]);
                  }
              }
              if (!qflag) {
                  printf("\n");
              }
          }
      } else {
          if (!qflag) {
              printf("Using q-form transform:\n");
          }
          for (i = 0; i < 4; i++) {
              for (j = 0; j < 4; j++) {
                  nii_xfm.m[i][j] = nii_ptr->qto_xyz.m[i][j];
                  if (!qflag) {
                      printf("%8.4f, ", nii_ptr->qto_xyz.m[i][j]);
                  }
              }
              if (!qflag) {
                  printf("\n");
              }
          }
      }

      /* Set up the correspondence between the file axes and the MINC
       * spatial axes. Each row contains the 'x', 'y', and 'z' components
       * along the right/left, anterior/posterior, or superior/inferior
       * axes (RAS). The "xspace" axis will be the one that has the largest
       * component in the RL direction, "yspace" refers to AP, and "zspace"
       * to SI. This tells us both how to convert the transform and how the
       * file data is arranged.
       */
      for ( i = 0; i < VIO_N_DIMENSIONS; i++) {
        int spatial_axis = VIO_X;
        char *dimname = MIxspace;
        float c_x = fabsf(nii_xfm.m[VIO_X][i]);
        float c_y = fabsf(nii_xfm.m[VIO_Y][i]);
        float c_z = fabsf(nii_xfm.m[VIO_Z][i]);
        if (c_y > c_x && c_y > c_z) {
          spatial_axis = VIO_Y;
          dimname = MIyspace;
        }
        if (c_z > c_x && c_z > c_y) {
          spatial_axis = VIO_Z;
          dimname = MIzspace;
        }
        spatial_axes[i] = spatial_axis;
        mnc_ordered_dim_names[i] = dimname;
        for (j = 0; j < 4; j++) {
          Transform_elem(mnc_xfm, j, spatial_axis) = nii_xfm.m[j][i];
        }
      }
      /* Put the final column in the right place.
       */
      for (i = 0; i < 4; i++) {
        Transform_elem(mnc_xfm, i, 3) = nii_xfm.m[i][3];
      }

      create_linear_transform(&mnc_linear_xfm, &mnc_xfm);
      convert_transform_to_starts_and_steps(&mnc_linear_xfm,
                                            VIO_N_DIMENSIONS,
                                            NULL,
                                            trivial_axes,
                                            dim_starts,
                                            dim_steps,
                                            dim_dircos);
    } else {
      /* No official transform was found (possibly this is an Analyze
       * file).  Just use some reasonable defaults.
       */
      mnc_ordered_dim_names[VIO_X] = MIxspace;
      mnc_ordered_dim_names[VIO_Y] = MIyspace;
      mnc_ordered_dim_names[VIO_Z] = MIzspace;
      dim_steps[VIO_X] = nii_ptr->dx;
      dim_steps[VIO_Y] = nii_ptr->dy;
      dim_steps[VIO_Z] = nii_ptr->dz;
      dim_starts[VIO_X] = -(nii_ptr->dx * nii_ptr->nx) / 2;
      dim_starts[VIO_Y] = -(nii_ptr->dy * nii_ptr->ny) / 2;
      dim_starts[VIO_Z] = -(nii_ptr->dz * nii_ptr->nz) / 2;

      /* Unlike the starts and steps, the direction cosines do NOT change
       * based upon the data orientation.
       */
      for (i = 0; i < VIO_N_DIMENSIONS; i++) {
          for (j = 0; j < VIO_N_DIMENSIONS; j++) {
              dim_dircos[i][j] = (i == j) ? 1.0 : 0.0;
          }
      }
    }

    /* Open the MINC file.  It should not already exist.
     */
    mnc_fd = micreate(out_str, NC_NOCLOBBER);
    if (mnc_fd < 0) {
        fprintf(stderr, "Can't create output file '%s'\n", out_str);
        return (-1);
    }

    /* Create the necessary dimensions in the minc file, starting from
     * the fastest-varying dimension.
     */

    img_ndims = 0;

    if (nii_ptr->nt > 1) {
        img_dimids[img_ndims] = ncdimdef(mnc_fd, MItime, nii_ptr->nt);
        img_count[img_ndims] = nii_ptr->nt;
        img_ndims++;

        r = micreate_std_variable(mnc_fd, MItime, NC_INT, 0, NULL);
        switch (nii_ptr->time_units) {
        case NIFTI_UNITS_UNKNOWN:
        case NIFTI_UNITS_SEC:
            time_step = nii_ptr->dt;
            time_start = nii_ptr->toffset;
            break;
        case NIFTI_UNITS_MSEC:
            time_step = nii_ptr->dt / 1000;
            time_start = nii_ptr->toffset / 1000;
            break;
        case NIFTI_UNITS_USEC:
            time_step = nii_ptr->dt / 1000000;
            time_start = nii_ptr->toffset / 1000000;
            break;
        default:
            fprintf(stderr, "Unknown time units value %d\n",
                    nii_ptr->time_units);
            break;
        }
        miattputdbl(mnc_fd, r, MIstart, time_start);
        miattputdbl(mnc_fd, r, MIstep, time_step);
        miattputstr(mnc_fd, r, MIunits, "s");
    }

    if (nii_ptr->nz > 1) {
        img_dimids[img_ndims] = ncdimdef(mnc_fd, mnc_ordered_dim_names[VIO_Z],
                                         nii_ptr->nz);
        img_count[img_ndims] = nii_ptr->nz;
        img_ndims++;

        r = micreate_std_variable(mnc_fd, mnc_ordered_dim_names[VIO_Z], NC_INT,
                                  0, NULL);
        miattputdbl(mnc_fd, r, MIstep, nii_ptr->dz);
        miattputstr(mnc_fd, r, MIunits, "mm");
    }

    if (nii_ptr->ny > 1) {
        img_dimids[img_ndims] = ncdimdef(mnc_fd, mnc_ordered_dim_names[VIO_Y],
                                         nii_ptr->ny);
        img_count[img_ndims] = nii_ptr->ny;
        img_ndims++;

        r = micreate_std_variable(mnc_fd, mnc_ordered_dim_names[VIO_Y], NC_INT,
                                  0, NULL);
        miattputdbl(mnc_fd, r, MIstep, nii_ptr->dy);
        miattputstr(mnc_fd, r, MIunits, "mm");
    }

    if (nii_ptr->nx > 1) {
        img_dimids[img_ndims] = ncdimdef(mnc_fd, mnc_ordered_dim_names[VIO_X],
                                         nii_ptr->nx);
        img_count[img_ndims] = nii_ptr->nx;
        img_ndims++;

        r = micreate_std_variable(mnc_fd, mnc_ordered_dim_names[VIO_X], NC_INT,
                                  0, NULL);
        miattputdbl(mnc_fd, r, MIstep, nii_ptr->dx);
        miattputstr(mnc_fd, r, MIunits, "mm");
    }

    if (nii_ptr->nu > 1) {
        img_dimids[img_ndims] = ncdimdef(mnc_fd, MIvector_dimension,
                                         nii_ptr->nu);
        img_count[img_ndims] = nii_ptr->nu;
        img_ndims++;
    }

    /* Create scalar image-min and image-max variables.
     */
    micreate_std_variable(mnc_fd, MIimagemax, NC_DOUBLE, 0, NULL);
    micreate_std_variable(mnc_fd, MIimagemin, NC_DOUBLE, 0, NULL);

    /* Create the group variables.
     */
    micreate_std_variable(mnc_fd, MIstudy, NC_INT, 0, NULL);
    if (strlen(nii_ptr->descrip) > 0 && strlen(nii_ptr->descrip) < 79 ) {
      int varid = micreate_std_variable(mnc_fd, MIpatient, NC_INT, 0, NULL);
      miattputstr(mnc_fd, varid, MIfull_name, nii_ptr->descrip);
    } else {
       micreate_std_variable(mnc_fd, MIpatient, NC_INT, 0, NULL);
    }
    micreate_std_variable(mnc_fd, MIacquisition, NC_INT, 0, NULL);

    /* Create the MINC image variable.  If we can't, there is no
     * further processing possible...
     */
    mnc_iid = micreate_std_variable(mnc_fd, MIimage, img_vtype, img_ndims,
                                    img_dimids);
    if (mnc_iid < 0) {
        fprintf(stderr, "Can't create the image variable\n");
        return (-1);
    }

    miattputstr(mnc_fd, mnc_iid, MIsigntype,
                (img_vsign) ? MI_SIGNED : MI_UNSIGNED);

    switch (nii_ptr->xyz_units) {
    case NIFTI_UNITS_METER:
        for (i = 0; i < VIO_N_DIMENSIONS; i++) {
            dim_starts[i] *= 1000;
            dim_steps[i] *= 1000;
        }
        break;
    case NIFTI_UNITS_MM:
    case NIFTI_UNITS_UNKNOWN:   /* Assume millimeters!! */
        break;
    case NIFTI_UNITS_MICRON:
        for (i = 0; i < VIO_N_DIMENSIONS; i++) {
            dim_starts[i] /= 1000;
            dim_steps[i] /= 1000;
        }
        break;
    default:
        fprintf(stderr, "Unknown XYZ units %d\n", nii_ptr->xyz_units);
        break;
    }

    if (!qflag) {
      printf("\n");
      printf("name   | start      | step    | cosines\n");
      for (i = VIO_N_DIMENSIONS-1; i >= 0; i--) {
        int k = spatial_axes[i];
        printf("%-6s | %10.4f | %7.4f | %7.4f %7.4f %7.4f\n",
               spatial_names[k],
               dim_starts[k],
               dim_steps[k],
               dim_dircos[k][VIO_X],
               dim_dircos[k][VIO_Y],
               dim_dircos[k][VIO_Z]);
      }
    }
    /* Now we write the spatial axis information to the file.  The starts,
     * steps, and cosines have to be associated with the correct spatial
     * axes.
     */
    for (i = 0; i < VIO_N_DIMENSIONS; i++) {
        int mnc_vid = ncvarid(mnc_fd, spatial_names[i]);

        miattputdbl(mnc_fd, mnc_vid, MIstart, dim_starts[i]);
        miattputdbl(mnc_fd, mnc_vid, MIstep, dim_steps[i]);
        ncattput(mnc_fd, mnc_vid, MIdirection_cosines, NC_DOUBLE,
                 VIO_N_DIMENSIONS, dim_dircos[i]);
    }

    /* Find the valid minimum and maximum of the data, in order to set the
     * global image minimum and image maximum properly.
     */
    if (rflag) {
        find_data_range(nii_ptr->datatype, nii_ptr->nvox, nii_ptr->data,
                        img_vrange);
    }

    if (nii_ptr->scl_slope != 0.0) {
        /* Convert slope/offset to min/max
         */
        img_rrange[0] = (img_vrange[0] * nii_ptr->scl_slope) + nii_ptr->scl_inter;
        img_rrange[1] = (img_vrange[1] * nii_ptr->scl_slope) + nii_ptr->scl_inter;
    }
    else {
        img_rrange[0] = img_vrange[0];
        img_rrange[1] = img_vrange[1];
    }

    ncattput(mnc_fd, mnc_iid, MIvalid_range, NC_DOUBLE, 2, img_vrange);
    miattputstr(mnc_fd, NC_GLOBAL, MIhistory, time_stamp(argc, argv));
    miattputstr(mnc_fd, mnc_iid, MIcomplete, MI_FALSE);

    /* Switch out of definition mode.
     */
    ncendef(mnc_fd);

    /* Finally, write the values of the image-min, image-max, and image
     * variables.
     */
    mivarput1(mnc_fd, ncvarid(mnc_fd, MIimagemin), img_start, NC_DOUBLE,
              MI_SIGNED, &img_rrange[0]);

    mivarput1(mnc_fd, ncvarid(mnc_fd, MIimagemax), img_start, NC_DOUBLE,
              MI_SIGNED, &img_rrange[1]);

    mivarput(mnc_fd, mnc_iid, img_start, img_count, img_mtype,
             (img_msign) ? MI_SIGNED : MI_UNSIGNED, nii_ptr->data);

    miattputstr(mnc_fd, mnc_iid, MIcomplete, MI_TRUE);
    miclose(mnc_fd);

    return (0);
}
