alpha_value = 0.5
text_size = 18
line_width = 1.5

calculatePosterior <- function(metrics, model) {
    # Determine which model component is the dead cell distribution,
    # since flexmix doesn't number the distributions deterministically.
    # Dead cell distribution should always have a higher y-intercept than
    # the alive cell distribution, no matter the kind of model used.
    intercept1 <- parameters(model, component = 1)[1]
    intercept2 <- parameters(model, component = 2)[1]
    if (intercept1 > intercept2) {
        dead_cell_dist <- 1
    } else {
        dead_cell_dist <- 2
    }
    
    # calculate posterior likelihood for each cell
    post <- posterior(model)
    metrics$probability_dead <- post[, dead_cell_dist]
    
    metrics
}

plotNoModel <- function(metrics, title = NULL, margin =T) {
    p <- ggplot(metrics, aes(x = detected, y = subsets_mito_percent)) +
        geom_point(colour = "#33adff", alpha = alpha_value) +
        labs(x = "Unique genes found", y = "Percent reads mitochondrial") +
        theme_bw()
    if (!is.null(title)){
        p <- p + ggtitle(title)
    }    
    if (margin) {
        p <- p + theme(text = element_text(size = text_size),
                       plot.margin = unit(c(1, 1, 1, 1), "cm"))
    }
    
    p
}

plotDistributions <- function(metrics, fitted_models, 
                              title = NULL, margin = T) {
    p <- ggplot(metrics, aes(x = detected, y = subsets_mito_percent)) +
        geom_point(colour = "#33adff", alpha = alpha_value) +
        geom_line(data = fitted_models, inherit.aes = F,
                  aes(x = V1, y = Comp.1), lwd = line_width) +
        geom_line(data = fitted_models, inherit.aes = F,
                  aes(x = V1, y = Comp.2), lwd = line_width) +
        labs(x = "Unique genes found", y = "Percent reads mitochondrial") +
        theme_bw() + ylim(0, max(metrics$subsets_mito_percent))
    if (!is.null(title)){
        p <- p + ggtitle(title)
    }    
    if (margin) {
        p <- p + theme(text = element_text(size = text_size),
                       plot.margin = unit(c(1, 1, 1, 1), "cm"))
    }
    
    p
}

plotPosterior <- function(metrics, fitted_models,
                          title = NULL, margin = T){
    p <- ggplot(metrics, aes(x = detected, y = subsets_mito_percent,
                              colour = probability_dead)) +
        labs(x = "Unique genes found", y = "Percent reads mitochondrial",
             colour = "Probability\ncompromised") +
        geom_point(alpha = alpha_value) +
        geom_line(data = fitted_models, inherit.aes = F,
                  aes(x = V1, y = Comp.1), lwd = line_width) +
        geom_line(data = fitted_models, inherit.aes = F, 
                  aes(x = V1, y = Comp.2), lwd = line_width) +
        theme_bw() + ylim(0, max(metrics$subsets_mito_percent))
    if (!is.null(title)){
        p <- p + ggtitle(title)
    }    
    if (margin) {
        p <- p + theme(text = element_text(size = text_size),
                       plot.margin = unit(c(1, 1, 1, 1), "cm"))
    }
    
    p
}

plotKeepOrToss <- function(metrics, title = NULL, margin = T) {
    p <- ggplot(metrics, aes(x = detected, y = subsets_mito_percent,
                             colour = keep)) +
        scale_color_manual(values = c("#999999", "#E69F00")) +
        labs(x = "Unique genes found", y = "Percent reads mitochondrial",
             colour = "Keep") +
        geom_point(alpha = alpha_value) + theme_bw()
    if (!is.null(title)){
        p <- p + ggtitle(title)
    }    
    if (margin) {
        p <- p + theme(text = element_text(size = text_size),
                       plot.margin = unit(c(1, 1, 1, 1), "cm"))
    }
    
    p
}