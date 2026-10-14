new_tensor <- function(x) {
  if(!is.array(x)) {
    stop("`x` must be an array.")
  }

  new_vctr(.data = x, dim = dim(x), class = "tensor")
}

tensor <- function(x, dim = NULL, draw_dim = NULL) {
  if (is.null(dim)) { #no input dims
    dim <- dim(x) #keep dims of input

    if (is.null(dim)) { #input was a vector
      dim <- length(x) #take whole length
    }
  }
  
  order <- length(dim)

  if(is.null(draw_dim)) {
    new_tensor(array(x, dim = c(1, dim)))
  }
  
  else { #place draw dim at front
    new_tensor(aperm(array(x, dim = dim), c(draw_dim, (1:order)[-draw_dim])))
  }
}

n_draws <- function(x) {
  vctrs::vec_size(x)
}

shape <- function(x) {
  dim(x)[-1L]
}

order <- function(x) {
  length(shape(x))
}

pull_draw <- function(x, i) {
  if (!inherits(x, "tensor")) {
    stop("`x` must be a tensor.", call. = FALSE)
  }

  i <- vctrs::vec_cast(i, integer())

  if (length(i) != 1L || is.na(i) || i < 1L || i > n_draws(x)) {
    stop("`i` must identify exactly one draw.", call. = FALSE)
  }

  sliced <- vctrs::vec_slice(x, i)

  array(as.vector(vctrs::vec_data(sliced)), dim = shape(x))
}

print.tensor <- function(
  x,
  ...,
  max_draws = getOption("tensortools.print_max_draws", 3L),
  max_slices = getOption("tensortools.print_max_slices", 4L)
) {
  size <- n_draws(x)
  shape <- shape(x)
  order <- order(x)
  type <- typeof(vctrs::vec_data(x))

  cat(
    "<tensor[", size, " draws]>\n",
    "Shape per draw: ", paste(shape, collapse = " \u00d7 "), "\n",
    "Type: ", type, "\n",
    sep = ""
  )

  draws_to_show <- seq_len(min(size, max_draws))

  for (i in draws_to_show) {
    cat("\nDraw ", i, ":\n", sep = "")

    draw <- pull_draw(x, i)

    if (order == 1L) {
      print(as.vector(draw), ...)
      next
    }

    if (
      order == 3L &&
      shape[3L] > max_slices
    ) {
      print(draw[, , seq_len(max_slices), drop = FALSE], ...)

      omitted <- shape[3L] - max_slices

      cat(
        "... ", omitted, " more ",
        if (omitted == 1L) "slice" else "slices",
        "\n",
        sep = ""
      )

      next
    }

    print(draw, ...)
  }

  if (size > max_draws) {
    omitted <- size - max_draws

    cat(
      "\n... ", omitted, " more ",
      if (omitted == 1L) "draw" else "draws",
      "\n",
      sep = ""
    )
  }

  invisible(x)
}
